!> This module exists to provide a wrapper around the MPI module. This is all because using the compiler flag -fdefault-integer-8 makes it harder to work with MPI.
module mpi_wrapper_module
   use mpi
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

   !> Subroutine wrapper for a non-blocking send and recieve operation
   subroutine isendrecv_mpi_wrapper(array_size, send_array, recv_array, array_start_index, count, &
      sendrecv_rank, comm, send_request_8, recv_request_8)

      integer, intent(in) :: array_size, array_start_index, count, sendrecv_rank, comm
      real, dimension(array_size), intent(inout) :: send_array, recv_array
      integer, intent(out) :: send_request_8, recv_request_8

      call MPI_ISEND(send_array(array_start_index), INT(count,kind=4), MPI_DOUBLE_PRECISION, INT(sendrecv_rank,kind=4),&
         neighbor_sendrecv_tag, INT(comm,kind=4), send_request_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_IRECV(recv_array(array_start_index), INT(count,kind=4), MPI_DOUBLE_PRECISION, INT(sendrecv_rank,kind=4),&
         neighbor_sendrecv_tag, INT(comm,kind=4), recv_request_4, ierr_4)
      call check_error_mpi(ierr_4)

      send_request_8 = send_request_4
      recv_request_8 = recv_request_4

   end subroutine isendrecv_mpi_wrapper

   !> Subroutine wrapper for a wait all operation
   subroutine waitall_mpi_wrapper(num_requests_8, requests_8)
      integer, intent(in) :: num_requests_8
      integer, dimension(num_requests_8), intent(inout) :: requests_8

      integer(kind=4), dimension(num_requests_8) :: requests_4

      requests_4 = requests_8

      call MPI_WAITALL(INT(num_requests_8,kind=4), requests_4, MPI_STATUSES_IGNORE, ierr_4)
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

      call MPI_CART_CREATE(INT(MPI_COMM_WORLD,kind=4), INT(ndims_8,kind=4), INT(dims_8,kind=4),&
         periods_4, reorder_4, comm_cart_4, ierr_4)
      call check_error_mpi(ierr_4)

      comm_cart_8 = comm_cart_4
   end subroutine create_cart_communicator_mpi_wrapper

   !> Get the Cartesian coordinates of the current process
   subroutine get_cart_coords_mpi_wrapper(comm_8, rank_8, ndims_8, coords_8)
      integer, intent(in) :: comm_8, rank_8, ndims_8

      integer(kind=4) :: coords_4(ndims_8)

      integer, intent(out) :: coords_8(ndims_8)

      call MPI_CART_COORDS(INT(comm_8,kind=4), INT(rank_8,kind=4), INT(ndims_8,kind=4), coords_4, ierr_4)
      call check_error_mpi(ierr_4)

      coords_8 = coords_4
   end subroutine get_cart_coords_mpi_wrapper

   !> Routine to get the rank of a process given its coordinates
   subroutine cart_rank_mpi_wrapper(comm_8, ndims_8, indices_8, rank_of_coords_8, error_8)
      integer, intent(in) :: comm_8, ndims_8
      integer, dimension(ndims_8), intent(in) :: indices_8

      integer, intent(out) :: rank_of_coords_8, error_8

      error_4 = 0
      call MPI_CART_RANK(INT(comm_8,kind=4), INT(indices_8,kind=4), rank_of_coords_4, ierr_4)
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
         call MPI_COMM_GET_ERRHANDLER(INT(comm_8,kind=4), original_errhandler_4, ierr_4)
         call check_error_mpi(ierr_4)
      endif

      ! Set the error handler for comm to return errors
      call MPI_COMM_SET_ERRHANDLER(INT(comm_8,kind=4), MPI_ERRORS_RETURN, ierr_4)
      call check_error_mpi(ierr_4)

   end subroutine change_MPI_COMM_errhandler_mpi_wrapper

   !> Restores the MPI_COMM error handler back to the original after finding neighbors
   subroutine original_MPI_COMM_errhandler_mpi_wrapper(comm_8)
      integer, intent(in) :: comm_8

      ! Check if the original error handler was saved
      if (original_errhandler_4 /= MPI_ERRHANDLER_NULL) then
         ! Restore the original error handler for comm_8
         call MPI_COMM_SET_ERRHANDLER(INT(comm_8,kind=4), original_errhandler_4, ierr_4)
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
   subroutine write_to_file_mpi_wrapper(filename, comm, buffer, count, datatype, offset)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: comm, count, datatype
      real, dimension(count), intent(in) :: buffer

      integer(kind=4) :: comm_4, count_4, datatype_4

      integer(kind=MPI_OFFSET_KIND), intent(in) :: offset
      integer(kind=4) :: file_handle_4
      integer(kind=4), dimension(MPI_STATUS_SIZE) :: status

      comm_4 = comm
      count_4 = count
      datatype_4 = datatype

      call MPI_FILE_OPEN(comm_4, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file_handle_4, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_FILE_WRITE_AT(file_handle_4, offset, buffer, count_4, datatype_4, status, ierr_4)
      call check_error_mpi(ierr_4)

      call MPI_FILE_CLOSE(file_handle_4, ierr_4)
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
