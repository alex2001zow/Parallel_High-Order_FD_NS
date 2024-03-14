!> This module exists to provide a wrapper around the MPI module. This is all because using the compiler flag -fdefault-integer-8 makes it harder to work with MPI.
module mpi_wrapper_module
   use mpi
   use utility_functions_module, only: sleeper_function
   implicit none

   private

   ! 4 byte integers because of OpenMPI
   integer(kind=4) :: rank_4, world_size_4, ierr_4, comm_cart_4, rank_of_coords_4, error_4
   integer(kind=4) :: original_errhandler_4 = MPI_ERRHANDLER_NULL

   integer(kind=4) :: send_request_4, recv_request_4
   integer(kind=4), parameter :: neighbor_sendrecv_tag = 99

   public :: initialize_mpi_wrapper, finalize_mpi_wrapper, create_cart_communicator_mpi_wrapper, &
      get_cart_coords_mpi_wrapper, change_MPI_COMM_errhandler_mpi_wrapper, original_MPI_COMM_errhandler_mpi_wrapper, &
      cart_rank_mpi_wrapper, isendrecv_mpi_wrapper, waitall_mpi_wrapper, free_communicator_mpi_wrapper, all_reduce_mpi_wrapper, &
      write_to_file_mpi_wrapper, check_openmpi_version

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

   subroutine isend_mpi_wrapper(array_size, send_array, array_start_index, count, sendrecv_rank, comm, send_request_8)

      integer, intent(in) :: array_size, array_start_index, count, sendrecv_rank, comm
      real, dimension(array_size), intent(inout) :: send_array
      integer, intent(out) :: send_request_8

      call MPI_ISEND(send_array(array_start_index), int(count,kind=4), MPI_DOUBLE_PRECISION, INT(sendrecv_rank,kind=4), &
         neighbor_sendrecv_tag, int(comm,kind=4), send_request_4, ierr_4)
      call check_error_mpi(ierr_4)

      send_request_8 = send_request_4

   end subroutine isend_mpi_wrapper

   subroutine irecv_mpi_wrapper(array_size, recv_array, array_start_index, count, sendrecv_rank, comm, recv_request_8)

      integer, intent(in) :: array_size, array_start_index, count, sendrecv_rank, comm
      real, dimension(array_size), intent(inout) :: recv_array
      integer, intent(out) :: recv_request_8

      call MPI_IRECV(recv_array(array_start_index), int(count,kind=4), MPI_DOUBLE_PRECISION, int(sendrecv_rank,kind=4), &
         neighbor_sendrecv_tag, int(comm,kind=4), recv_request_4, ierr_4)
      call check_error_mpi(ierr_4)

      recv_request_8 = recv_request_4

   end subroutine irecv_mpi_wrapper

   !> Subroutine wrapper for a non-blocking send and recieve operations
   subroutine isendrecv_mpi_wrapper(array_size, send_array, recv_array, array_start_index, count, &
      sendrecv_rank, comm, send_request_8, recv_request_8)

      integer, intent(in) :: array_size, array_start_index, count, sendrecv_rank, comm
      real, dimension(array_size), intent(inout) :: send_array, recv_array
      integer, intent(out) :: send_request_8, recv_request_8

      ! Declaration section remains the same
      call isend_mpi_wrapper(array_size, send_array, array_start_index, count, sendrecv_rank, comm, send_request_8)
      call irecv_mpi_wrapper(array_size, recv_array, array_start_index, count, sendrecv_rank, comm, recv_request_8)

   end subroutine isendrecv_mpi_wrapper

   !> Subroutine wrapper for a wait all operation
   subroutine waitall_mpi_wrapper(num_requests_8, requests_8)
      integer, intent(in) :: num_requests_8
      integer, dimension(num_requests_8), intent(inout) :: requests_8

      integer(kind=4), dimension(num_requests_8) :: requests_4

      requests_4 = requests_8

      call MPI_WAITALL(int(num_requests_8,kind=4), requests_4, MPI_STATUSES_IGNORE, ierr_4)
      call check_error_mpi(ierr_4)

      requests_8 = requests_4
   end subroutine waitall_mpi_wrapper

   !> This subroutine creates a cartesian communicator
   subroutine create_cart_communicator_mpi_wrapper(ndims_8, dims_8, comm_cart_8)
      integer, intent(in) :: ndims_8
      integer, dimension(ndims_8), intent(in) :: dims_8

      logical(kind=4) :: periods_4(ndims_8)
      logical(kind=4) :: reorder_4 = .TRUE.

      integer, intent(out) :: comm_cart_8

      integer :: ii

      do ii = 1, ndims_8
         periods_4(ii) = .FALSE.
      end do

      call MPI_CART_CREATE(INT(MPI_COMM_WORLD,kind=4), int(ndims_8,kind=4), int(dims_8,kind=4),&
         periods_4, reorder_4, comm_cart_4, ierr_4)
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
      integer, dimension(ndims_8), intent(in) :: indices_8

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

   !> Write the block to the file using MPI-IO
   subroutine write_to_file_mpi_wrapper(ndims, solution_filename, comm, num_blocks, block_length, block_stride, &
      starting_offset, extent, global_dims, block_dims, block_begin, solution_buffer, solution_offset, solution_count)

      integer, intent(in) :: ndims, num_blocks, block_length, block_stride, starting_offset, extent
      character(len=*), intent(in) :: solution_filename
      integer, intent(in) :: comm, solution_count
      integer, dimension(ndims), intent(in) :: global_dims, block_dims, block_begin
      real, dimension(solution_count), intent(in) :: solution_buffer
      integer(kind=MPI_OFFSET_KIND), intent(in) :: solution_offset


      integer(kind=4) :: ndims_4, num_blocks_4, block_length_4, block_stride_4
      integer(kind=MPI_ADDRESS_KIND) :: starting_offset_4, extent_4
      integer(kind=4) :: comm_4, solution_count_4
      integer(kind=4), dimension(ndims) :: global_dims_4, block_dims_4, block_begin_4
      integer(kind=4) :: block_type, reshaped_block_type, file_type
      integer(kind=4) :: solution_fh_4
      integer(kind=4), dimension(MPI_STATUS_SIZE) :: solution_status

      ndims_4 = ndims
      num_blocks_4 = num_blocks
      block_length_4 = block_length
      block_stride_4 = block_stride
      starting_offset_4 = starting_offset
      extent_4 = extent

      comm_4 = comm
      solution_count_4 = 1

      global_dims_4 = global_dims
      block_dims_4 = block_dims
      block_begin_4 = block_begin

      ! Specify the block setup as vector with stride
      call MPI_TYPE_VECTOR(num_blocks_4, block_length_4, block_stride_4, MPI_DOUBLE_PRECISION, block_type, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_TYPE_COMMIT(block_type, ierr_4)
      call check_error_mpi(ierr_4)

      ! Get extent
      call MPI_TYPE_GET_EXTENT(block_type, starting_offset_4, extent_4, ierr_4)

      ! Reshape the block to start and end at the correct places
      call MPI_TYPE_CREATE_RESIZED(block_type, starting_offset_4, extent_4, reshaped_block_type, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_TYPE_COMMIT(reshaped_block_type, ierr_4)
      call check_error_mpi(ierr_4)

      ! Specify the block dimensions in the global space.
      call MPI_TYPE_CREATE_SUBARRAY(ndims_4, global_dims_4, block_dims_4, block_begin_4, &
         MPI_ORDER_C, MPI_DOUBLE_PRECISION, file_type, ierr_4) ! The order is C because the matrix has been indexed using a global index with row major order. Should not effect performance.
      call check_error_mpi(ierr_4)

      call MPI_TYPE_COMMIT(file_type, ierr_4)
      call check_error_mpi(ierr_4)

      ! Write the solution to the file
      call MPI_FILE_OPEN(comm_4, solution_filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &
         solution_fh_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_FILE_SET_VIEW(solution_fh_4, solution_offset, MPI_DOUBLE_PRECISION, file_type, "native", MPI_INFO_NULL, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_FILE_WRITE_ALL(solution_fh_4, solution_buffer, solution_count_4, reshaped_block_type, solution_status, ierr_4)
      call check_error_mpi(ierr_4)

      ! Clean up
      call MPI_FILE_CLOSE(solution_fh_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_TYPE_FREE(reshaped_block_type, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_TYPE_FREE(block_type, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_TYPE_FREE(file_type, ierr_4)
      call check_error_mpi(ierr_4)
   end subroutine write_to_file_mpi_wrapper

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
