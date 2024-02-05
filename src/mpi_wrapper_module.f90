!> This module exists to provide a wrapper around the MPI module. This is all because using the compiler flag -fdefault-integer-8 makes it harder to work with MPI.
module mpi_wrapper_module
   use mpi
   implicit none

   private
   !character(len=MPI_MAX_ERROR_STRING) :: error_string
   !integer(kind=4) :: ierr_string, error_len

   ! 4 byte integers because of OpenMPI
   integer(kind=4) :: rank_4, world_size_4, ierr_4, comm_cart_4, rank_of_coords_4, error_4
   integer(kind=4) :: original_errhandler_4 = MPI_ERRHANDLER_NULL

   integer(kind=4) :: send_request_4, recv_request_4
   integer(kind=4), parameter :: neighbor_sendrecv_tag = 99

   public :: initialize_mpi_wrapper, finalize_mpi_wrapper, create_cart_communicator_mpi_wrapper, &
      get_cart_coords_mpi_wrapper, change_MPI_COMM_errhandler_mpi_wrapper, original_MPI_COMM_errhandler_mpi_wrapper, &
      cart_rank_mpi_wrapper, isendrecv_mpi_wrapper

contains

   !> This subroutine initializes MPI and gets the rank and world size
   subroutine initialize_mpi_wrapper(rank_8, world_size_8, ierr_8)
      integer, intent(out) :: rank_8, world_size_8, ierr_8

      call MPI_INIT(ierr_4)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank_4, ierr_4)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size_4, ierr_4)

      rank_8 = rank_4
      world_size_8 = world_size_4
      ierr_8 = ierr_4
   end subroutine initialize_mpi_wrapper

   !> This subroutine finalizes MPI
   subroutine finalize_mpi_wrapper(ierr_8)
      integer, intent(out) :: ierr_8

      call MPI_FINALIZE(ierr_4)

      ierr_8 = ierr_4
   end subroutine finalize_mpi_wrapper

   !> Subroutine wrapper for a non-blocking send operation
   subroutine isendrecv_mpi_wrapper(array_size, send_array, recv_array, array_start_index, count, &
      sendrecv_rank, comm, send_request_8, recv_request_8, ierr_8)

      integer, intent(in) :: array_size, array_start_index, count, sendrecv_rank, comm
      real, dimension(array_size), intent(inout) :: send_array, recv_array
      integer, intent(out) :: send_request_8, recv_request_8, ierr_8

      call MPI_ISEND(send_array(array_start_index), INT(count,kind=4), MPI_DOUBLE_PRECISION, INT(sendrecv_rank,kind=4),&
         neighbor_sendrecv_tag, INT(comm,kind=4), send_request_4, ierr_4)

      call MPI_IRECV(recv_array(array_start_index), INT(count,kind=4), MPI_DOUBLE_PRECISION, INT(sendrecv_rank,kind=4),&
         neighbor_sendrecv_tag, INT(comm,kind=4), recv_request_4, ierr_4)

      send_request_8 = send_request_4
      recv_request_8 = recv_request_4
      ierr_8 = ierr_4

   end subroutine isendrecv_mpi_wrapper

   !> This subroutine creates a cartesian communicator
   subroutine create_cart_communicator_mpi_wrapper(ndims_8, dims_8, comm_cart_8, ierr_8)
      integer, intent(in) :: ndims_8
      integer, dimension(ndims_8), intent(in) :: dims_8

      logical(kind=4) :: periods_4(ndims_8)
      logical(kind=4) :: reorder_4 = .TRUE.

      integer, intent(out) :: comm_cart_8, ierr_8

      integer :: ii

      do ii = 1, ndims_8
         periods_4(ii) = .FALSE.
      end do

      call MPI_Cart_create(INT(MPI_COMM_WORLD,kind=4), INT(ndims_8,kind=4), INT(dims_8,kind=4),&
         periods_4, reorder_4, comm_cart_4, ierr_4)

      comm_cart_8 = comm_cart_4
      ierr_8 = ierr_4

   end subroutine create_cart_communicator_mpi_wrapper

   !> Get the Cartesian coordinates of the current process
   subroutine get_cart_coords_mpi_wrapper(comm_8, rank_8, ndims_8, coords_8, ierr_8)
      integer, intent(in) :: comm_8, rank_8, ndims_8

      integer(kind=4) :: coords_4(ndims_8)

      integer, intent(out) :: coords_8(ndims_8), ierr_8

      call MPI_Cart_coords(INT(comm_8,kind=4), INT(rank_8,kind=4), INT(ndims_8,kind=4), coords_4, ierr_4)

      coords_8 = coords_4
      ierr_8 = ierr_4
   end subroutine get_cart_coords_mpi_wrapper

   !> Routine to get the rank of a process given its coordinates
   subroutine cart_rank_mpi_wrapper(comm_8, ndims_8, indices_8, rank_of_coords_8, error_8, ierr_8)
      integer, intent(in) :: comm_8, ndims_8
      integer, dimension(ndims_8), intent(in) :: indices_8

      integer, intent(out) :: rank_of_coords_8, error_8, ierr_8

      error_4 = 0
      call MPI_Cart_rank(INT(comm_8,kind=4), INT(indices_8,kind=4), rank_of_coords_4, ierr_4)
      if(ierr_4 /= MPI_SUCCESS) then
         error_4 = 1
      endif

      rank_of_coords_8 = rank_of_coords_4
      error_8 = error_4
      ierr_8 = ierr_4
   end subroutine cart_rank_mpi_wrapper

   !> Routine to change the MPI_COMM error handler so we can find neighbors without crashing. We restore the original using original_MPI_COMM_errhandler() when done
   subroutine change_MPI_COMM_errhandler_mpi_wrapper(comm_8, ierr_8)
      integer, intent(in) :: comm_8

      integer, intent(out) :: ierr_8

      ! Check if the original error handler has already been saved
      if (original_errhandler_4 == MPI_ERRHANDLER_NULL) then
         ! Get the current error handler for comm
         call MPI_Comm_get_errhandler(INT(comm_8,kind=4), original_errhandler_4, ierr_4)
         if (ierr_4 /= MPI_SUCCESS) then
            print *, "Error getting original MPI error handler."
            return
         endif
      endif

      ! Set the error handler for comm to return errors
      call MPI_Comm_set_errhandler(INT(comm_8,kind=4), MPI_ERRORS_RETURN, ierr_4)
      if (ierr_4 /= MPI_SUCCESS) then
         print *, "Error setting MPI error handler to MPI_ERRORS_RETURN."
      endif

      ierr_8 = ierr_4
   end subroutine change_MPI_COMM_errhandler_mpi_wrapper

   !> Restores the MPI_COMM error handler back to the original after finding neighbors
   subroutine original_MPI_COMM_errhandler_mpi_wrapper(comm_8, ierr_8)
      integer, intent(in) :: comm_8

      integer, intent(out) :: ierr_8

      ! Check if the original error handler was saved
      if (original_errhandler_4 /= MPI_ERRHANDLER_NULL) then
         ! Restore the original error handler for parameters%cart_comm
         call MPI_Comm_set_errhandler(INT(comm_8,kind=4), original_errhandler_4, ierr_4)
         if (ierr_4 /= MPI_SUCCESS) then
            print *, "Error restoring original MPI error handler."
         endif
         ! Reset the original_errhandler variable
         original_errhandler_4 = MPI_ERRHANDLER_NULL
      else
         print *, "Original MPI error handler not saved."
      endif

      ierr_8 = ierr_4
   end subroutine original_MPI_COMM_errhandler_mpi_wrapper

end module mpi_wrapper_module
