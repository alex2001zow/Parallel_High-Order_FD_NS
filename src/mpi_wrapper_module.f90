!> This module exists to provide a wrapper around the MPI module. This is all because using the compiler flag -fdefault-integer-8 makes it harder to work with MPI.
module mpi_wrapper_module
   use mpi
   use utility_functions_module, only: sleeper_function
   implicit none

   private

   type block_data_layout_type
      integer :: number_of_existing_neighbors
      integer, dimension(:), allocatable :: neighbor_ranks
      integer(kind=4) :: elements_per_index_type_4, global_block_type_4, true_block_type_4
      integer(kind=MPI_OFFSET_KIND) :: writing_offset
      integer(kind=4), dimension(:), allocatable :: send_type_4, recv_type_4
      integer(kind=4), dimension(:), allocatable :: send_tag_4, recv_tag_4
   end type block_data_layout_type

   ! 4 byte integers because of OpenMPI
   integer(kind=4) :: rank_4, world_size_4, ierr_4, comm_cart_4, rank_of_coords_4, error_4
   integer(kind=4) :: original_errhandler_4 = MPI_ERRHANDLER_NULL

   integer(kind=4) :: send_request_4, recv_request_4
   integer(kind=4), parameter :: neighbor_sendrecv_tag = 99

   public :: block_data_layout_type, define_block_layout, create_send_recv_layout, sendrecv_data_to_neighbors, &
      free_block_layout_type, write_block_data_to_file
   public :: initialize_mpi_wrapper, finalize_mpi_wrapper, create_cart_communicator_mpi_wrapper, &
      get_cart_coords_mpi_wrapper, change_MPI_COMM_errhandler_mpi_wrapper, original_MPI_COMM_errhandler_mpi_wrapper, &
      cart_rank_mpi_wrapper, waitall_mpi_wrapper, free_communicator_mpi_wrapper, all_reduce_mpi_wrapper, &
      check_openmpi_version

contains

   !> This subroutine initializes MPI and gets the rank and world size
   subroutine initialize_mpi_wrapper(rank_8, world_size_8)
      integer, intent(out) :: rank_8, world_size_8

      call MPI_INIT(ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_COMM_RANK(MPI_COMM_WORLD, rank_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size_4, ierr_4)
      call check_error_mpi(ierr_4)

      rank_8 = rank_4
      world_size_8 = world_size_4
   end subroutine initialize_mpi_wrapper

   !> This subroutine finalizes MPI
   subroutine finalize_mpi_wrapper()

      call MPI_FINALIZE(ierr_4)
      call check_error_mpi(ierr_4)

   end subroutine finalize_mpi_wrapper

   !> Subroutine wrapper for a wait all operation
   subroutine waitall_mpi_wrapper(num_requests_8, requests_8)
      integer, intent(in) :: num_requests_8
      integer, dimension(:), intent(inout) :: requests_8

      integer(kind=4), dimension(num_requests_8) :: requests_4

      requests_4 = requests_8

      call MPI_WAITALL(int(num_requests_8,kind=4), requests_4, MPI_STATUSES_IGNORE, ierr_4)
      call check_error_mpi(ierr_4)

      requests_8 = requests_4
   end subroutine waitall_mpi_wrapper

   !> This subroutine creates a cartesian communicator. We need to incoropate the periods as an input.
   subroutine create_cart_communicator_mpi_wrapper(ndims_8, dims_8, periods_8, reorder_8, comm_cart_8)
      integer, intent(in) :: ndims_8
      logical, intent(in) :: reorder_8
      integer, dimension(:), intent(in) :: dims_8
      logical, dimension(:), intent(in) :: periods_8
      integer, intent(out) :: comm_cart_8

      call MPI_CART_CREATE(int(MPI_COMM_WORLD,kind=4), int(ndims_8,kind=4), int(dims_8,kind=4),&
         logical(periods_8,kind=4), logical(reorder_8,kind=4), comm_cart_4, ierr_4)
      call check_error_mpi(ierr_4)

      comm_cart_8 = comm_cart_4
   end subroutine create_cart_communicator_mpi_wrapper

   !> Get the Cartesian coordinates of the current process
   subroutine get_cart_coords_mpi_wrapper(comm_8, rank_8, ndims_8, coords_8)
      integer, intent(in) :: comm_8, rank_8, ndims_8

      integer(kind=4) :: coords_4(ndims_8)

      integer, intent(out) :: coords_8(ndims_8)

      call MPI_CART_COORDS(int(comm_8,kind=4), int(rank_8,kind=4), int(ndims_8,kind=4), coords_4, ierr_4)
      call check_error_mpi(ierr_4)

      coords_8 = coords_4
   end subroutine get_cart_coords_mpi_wrapper

   !> Routine to get the rank of a process given its coordinates
   subroutine cart_rank_mpi_wrapper(comm_8, ndims_8, indices_8, rank_of_coords_8, error_8)
      integer, intent(in) :: comm_8, ndims_8
      integer, dimension(:), intent(in) :: indices_8

      integer, intent(out) :: rank_of_coords_8, error_8

      error_4 = 0
      call MPI_CART_RANK(int(comm_8,kind=4), int(indices_8,kind=4), rank_of_coords_4, ierr_4)
      !call check_error_mpi(ierr_4) ! We do not want to check for errors here because we want to return an error code to the caller
      if(ierr_4 /= MPI_SUCCESS) then
         error_4 = 1
      endif

      rank_of_coords_8 = rank_of_coords_4
      error_8 = error_4
   end subroutine cart_rank_mpi_wrapper

   !> Free a Cartesian communicator
   subroutine free_communicator_mpi_wrapper(comm_8)
      integer, intent(in) :: comm_8
      integer(kind=4) :: comm_4

      comm_4 = comm_8

      call MPI_COMM_FREE(comm_4, ierr_4)
      call check_error_mpi(ierr_4)
   end subroutine free_communicator_mpi_wrapper

   !> Routine to change the MPI_COMM error handler so we can find neighbors without crashing. We restore the original using original_MPI_COMM_errhandler() when done
   subroutine change_MPI_COMM_errhandler_mpi_wrapper(comm_8)
      integer, intent(in) :: comm_8

      ! Check if the original error handler has already been saved
      if (original_errhandler_4 == MPI_ERRHANDLER_NULL) then
         ! Get the current error handler for comm
         call MPI_COMM_GET_ERRHANDLER(int(comm_8,kind=4), original_errhandler_4, ierr_4)
         call check_error_mpi(ierr_4)
      endif

      ! Set the error handler for comm to return errors
      call MPI_COMM_SET_ERRHANDLER(int(comm_8,kind=4), MPI_ERRORS_RETURN, ierr_4)
      call check_error_mpi(ierr_4)

   end subroutine change_MPI_COMM_errhandler_mpi_wrapper

   !> Restores the MPI_COMM error handler back to the original after finding neighbors
   subroutine original_MPI_COMM_errhandler_mpi_wrapper(comm_8)
      integer, intent(in) :: comm_8

      ! Check if the original error handler was saved
      if (original_errhandler_4 /= MPI_ERRHANDLER_NULL) then
         ! Restore the original error handler for comm_8
         call MPI_COMM_SET_ERRHANDLER(int(comm_8,kind=4), original_errhandler_4, ierr_4)
         call check_error_mpi(ierr_4)

         ! Reset the original_errhandler variable
         original_errhandler_4 = MPI_ERRHANDLER_NULL
      else
         print *, "Original MPI error handler not saved."
      endif

   end subroutine original_MPI_COMM_errhandler_mpi_wrapper

   !> This subroutine is a wrapper for the MPI_ALLREDUCE function
   subroutine all_reduce_mpi_wrapper(sendbuf, recvbuf, count, datatype, op, comm)
      integer, intent(in) :: count, datatype, op, comm
      real, dimension(count), intent(in) :: sendbuf
      real, dimension(count), intent(out) :: recvbuf

      integer(kind=4) :: count_4, datatype_4, op_4, comm_4

      count_4 = count
      datatype_4 = datatype
      op_4 = op
      comm_4 = comm

      call MPI_ALLREDUCE(sendbuf, recvbuf, count_4, datatype_4, op_4, comm_4, ierr_4)
      call check_error_mpi(ierr_4)
   end subroutine all_reduce_mpi_wrapper

   !> Subroutine to define the whole block, the true block and the offset for writing to a file.
   subroutine define_block_layout(ndims, elements_per_index, grid_size, &
      ghost_begin_c, ghost_end_c, &
      global_begin_c, global_dims, &
      extended_global_begin_c, extended_global_dims, &
      block_begin_c, block_dims, &
      extended_block_begin_c, extended_block_dims, &
      block_data_layout)
      integer :: ndims, elements_per_index
      integer, dimension(ndims), intent(in) :: grid_size, ghost_begin_c, ghost_end_c
      integer, dimension(ndims), intent(in) :: global_begin_c, global_dims, extended_global_begin_c, extended_global_dims
      integer, dimension(ndims), intent(in) :: block_begin_c, block_dims, extended_block_begin_c, extended_block_dims
      type(block_data_layout_type), intent(inout) :: block_data_layout

      integer(kind=4) :: ndims_4, elements_per_index_4
      integer(kind=4), dimension(ndims) :: fullsizes_4, subsizes_4, starts_4

      ndims_4 = ndims
      elements_per_index_4 = elements_per_index

      ! Define the number of elements for each index in the block
      call MPI_TYPE_CONTIGUOUS(elements_per_index_4, MPI_REAL8, block_data_layout%elements_per_index_type_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_TYPE_COMMIT(block_data_layout%elements_per_index_type_4, ierr_4)
      call check_error_mpi(ierr_4)

      ! Define the full block with stencil points and ghost points In the global block space.
      fullsizes_4 = grid_size
      subsizes_4 = global_dims
      starts_4 = global_begin_c
      where(starts_4 /= 0)
         starts_4 = starts_4 - ghost_begin_c !- 1 ! Change this depending on if we have ghost points or not and if we are blocking
      end where
      !print *, "Global block: ", "fullsizes_4: ", fullsizes_4, "subsizes_4: ", subsizes_4, "starts_4: ", starts_4
      call MPI_TYPE_CREATE_SUBARRAY(ndims_4, fullsizes_4, subsizes_4, starts_4, MPI_ORDER_FORTRAN, &
         block_data_layout%elements_per_index_type_4, block_data_layout%global_block_type_4, ierr_4)

      call MPI_TYPE_COMMIT(block_data_layout%global_block_type_4, ierr_4)
      call check_error_mpi(ierr_4)

      ! Define the true block without stencil points and ghost points. This is the real block that we want to write to a file.
      fullsizes_4 = extended_block_dims
      subsizes_4 = block_dims
      starts_4 = block_begin_c
      !print *, "Local block: ", "fullsizes_4: ", fullsizes_4, "subsizes_4: ", subsizes_4, "starts_4: ", starts_4
      call MPI_TYPE_CREATE_SUBARRAY(ndims_4, fullsizes_4, subsizes_4, starts_4, &
         MPI_ORDER_FORTRAN, block_data_layout%elements_per_index_type_4, block_data_layout%true_block_type_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_TYPE_COMMIT(block_data_layout%true_block_type_4, ierr_4)
      call check_error_mpi(ierr_4)

      block_data_layout%writing_offset = 0

   end subroutine define_block_layout

   !> Subroutine to define the neighbor send and recieve types
   subroutine create_send_recv_layout(ndims, neighbor_index, full_block_dims, &
      begin_send, end_send, begin_recv, end_recv, block_data_layout)
      integer :: ndims, neighbor_index
      integer, dimension(ndims), intent(in) :: full_block_dims
      integer, dimension(ndims), intent(in) :: begin_send, end_send, begin_recv, end_recv
      type(block_data_layout_type), intent(inout) :: block_data_layout

      integer(kind=4) :: ndims_4
      integer(kind=4), dimension(ndims) :: fullsizes_4, subsizes_4, starts_4

      ndims_4 = ndims

      fullsizes_4 = full_block_dims
      subsizes_4 = end_send - begin_send
      starts_4 = begin_send
      !print *, "Send array: ", "fullsizes_4: ", fullsizes_4, "subsizes_4: ", subsizes_4, "starts_4: ", starts_4
      call MPI_TYPE_CREATE_SUBARRAY(ndims_4, fullsizes_4, subsizes_4, starts_4, MPI_ORDER_FORTRAN, &
         block_data_layout%elements_per_index_type_4, block_data_layout%send_type_4(neighbor_index), ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_TYPE_COMMIT(block_data_layout%send_type_4(neighbor_index), ierr_4)
      call check_error_mpi(ierr_4)

      fullsizes_4 = full_block_dims
      subsizes_4 = end_recv - begin_recv
      starts_4 = begin_recv
      !print *, "Recv array: ", "fullsizes_4: ", fullsizes_4, "subsizes_4: ", subsizes_4, "starts_4: ", starts_4
      call MPI_TYPE_CREATE_SUBARRAY(ndims_4, fullsizes_4, subsizes_4, starts_4, MPI_ORDER_FORTRAN, &
         block_data_layout%elements_per_index_type_4, block_data_layout%recv_type_4(neighbor_index), ierr_4)

      call MPI_TYPE_COMMIT(block_data_layout%recv_type_4(neighbor_index), ierr_4)
      call check_error_mpi(ierr_4)

   end subroutine create_send_recv_layout

   !> Subroutine to send and recieve data to and from neighbors using the send and recieve types. Split this into two subroutines to make it easier to use. And to open the receive then compute and then send.
   subroutine sendrecv_data_to_neighbors(block_data_layout, matrix, comm_cart, &
      neighbor_rank, neighbor_index, send_request, recv_request)
      type(block_data_layout_type), intent(in) :: block_data_layout
      real, dimension(:), intent(inout) :: matrix
      integer, intent(in) :: comm_cart, neighbor_rank, neighbor_index
      integer, intent(out) :: send_request, recv_request

      integer(kind=4) :: neighbor_rank_4, neighbor_index_4

      comm_cart_4 = comm_cart
      neighbor_rank_4 = neighbor_rank
      neighbor_index_4 = neighbor_index

      call MPI_ISEND(matrix, int(1,kind=4), block_data_layout%send_type_4(neighbor_index_4), neighbor_rank_4, &
         block_data_layout%send_tag_4(neighbor_index_4), comm_cart_4, send_request_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_IRECV(matrix, int(1,kind=4), block_data_layout%recv_type_4(neighbor_index_4), neighbor_rank_4, &
         block_data_layout%recv_tag_4(neighbor_index_4), comm_cart_4, recv_request_4, ierr_4)
      call check_error_mpi(ierr_4)

      send_request = send_request_4
      recv_request = recv_request_4

   end subroutine sendrecv_data_to_neighbors

   !> Subroutine to free the block layout types
   subroutine free_block_layout_type(block_data_layout)
      type(block_data_layout_type), intent(inout) :: block_data_layout

      integer :: global_index

      if (allocated(block_data_layout%neighbor_ranks)) deallocate(block_data_layout%neighbor_ranks)

      call MPI_TYPE_FREE(block_data_layout%elements_per_index_type_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_TYPE_FREE(block_data_layout%global_block_type_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_TYPE_FREE(block_data_layout%true_block_type_4, ierr_4)
      call check_error_mpi(ierr_4)

      ! Free the neighbor send and recieve types.
      do global_index = 1, block_data_layout%number_of_existing_neighbors
         call MPI_TYPE_FREE(block_data_layout%send_type_4(global_index), ierr_4)
         call check_error_mpi(ierr_4)

         call MPI_TYPE_FREE(block_data_layout%recv_type_4(global_index), ierr_4)
         call check_error_mpi(ierr_4)

      end do

      if (allocated(block_data_layout%send_type_4)) deallocate(block_data_layout%send_type_4)
      if (allocated(block_data_layout%recv_type_4)) deallocate(block_data_layout%recv_type_4)
      if (allocated(block_data_layout%send_tag_4)) deallocate(block_data_layout%send_tag_4)
      if (allocated(block_data_layout%recv_tag_4)) deallocate(block_data_layout%recv_tag_4)
   end subroutine free_block_layout_type

   ! Subroutine to write block data to a file
   subroutine write_block_data_to_file(block_data_layout, solution_filename, comm, block_data)
      type(block_data_layout_type), intent(in) :: block_data_layout
      character(len=*), intent(in) :: solution_filename
      integer, intent(in) :: comm
      real, dimension(:), intent(in) :: block_data

      integer(kind=4) :: comm_4, solution_fh_4, elements_to_write_4
      integer(kind=MPI_OFFSET_KIND) :: offset_4
      integer(kind=4), dimension(MPI_STATUS_SIZE) :: solution_status

      comm_4 = comm
      offset_4 = block_data_layout%writing_offset ! Is zero
      elements_to_write_4 = 1

      ! Write the solution to the file
      call MPI_FILE_OPEN(comm_4, solution_filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &
         solution_fh_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_FILE_SET_VIEW(solution_fh_4, offset_4, MPI_DOUBLE_PRECISION, &
         block_data_layout%global_block_type_4, "native", MPI_INFO_NULL, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_FILE_WRITE_ALL(solution_fh_4, real(block_data,kind=8), elements_to_write_4, &
         block_data_layout%true_block_type_4, solution_status, ierr_4)
      call check_error_mpi(ierr_4)

      ! Clean up
      call MPI_FILE_CLOSE(solution_fh_4, ierr_4)
      call check_error_mpi(ierr_4)

   end subroutine write_block_data_to_file

   !> This subroutine checks for an MPI error and aborts if there is one
   subroutine check_error_mpi(ierr_4_input)
      integer(kind=4), intent(in) :: ierr_4_input
      character(len=MPI_MAX_ERROR_STRING) :: error_string
      integer(kind=4) :: error_len, ierr_4_temp

      if (ierr_4_input /= MPI_SUCCESS) then
         call MPI_ERROR_STRING(ierr_4_input, error_string, error_len, ierr_4_temp)
         print *, 'MPI Error: ', trim(error_string)
         call MPI_ABORT(MPI_COMM_WORLD, ierr_4_input, ierr_4_temp)
      end if
   end subroutine check_error_mpi

   !> Subroutine to check OPENMPI version
   subroutine check_openmpi_version()
      integer(kind=4) :: major_version, minor_version
      character(len=100) :: version_string

      call MPI_GET_VERSION(major_version, minor_version, ierr_4)
      call check_error_mpi(ierr_4)

      write(version_string, '(I3,".",I3,".",I3)') major_version, minor_version
      print *, 'OpenMPI version: ', trim(version_string)
   end subroutine check_openmpi_version

end module mpi_wrapper_module
