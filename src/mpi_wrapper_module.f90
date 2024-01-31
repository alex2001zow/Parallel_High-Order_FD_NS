!> This module exists to provide a wrapper around the MPI module. This is all because using the compiler flag -fdefault-integer-8 makes it harder to work with MPI.
module mpi_wrapper_module
   use mpi
   implicit none

   private
   !character(len=MPI_MAX_ERROR_STRING) :: error_string
   !integer(kind=4) :: ierr_string, error_len

   integer(kind=4) :: original_errhandler = MPI_ERRHANDLER_NULL

   public :: initialize_mpi_wrapper, finalize_mpi_wrapper, create_cart_communicator_mpi_wrapper, &
      get_cart_coords_mpi_wrapper, change_MPI_COMM_errhandler_mpi_wrapper, original_MPI_COMM_errhandler_mpi_wrapper, &
      mpi_cart_rank_mpi_wrapper

contains

   !> This subroutine initializes MPI and gets the rank and world size
   subroutine initialize_mpi_wrapper(rank, world_size, ierr)
      integer(kind=4), intent(out) :: rank, world_size, ierr

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, ierr)
   end subroutine initialize_mpi_wrapper

   !> This subroutine finalizes MPI
   subroutine finalize_mpi_wrapper(ierr)
      integer(kind=4), intent(out) :: ierr

      call MPI_FINALIZE(ierr)
   end subroutine finalize_mpi_wrapper

   !> This subroutine creates a cartesian communicator
   subroutine create_cart_communicator_mpi_wrapper(ndims, dims, periods, comm_cart, ierr)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims
      logical, dimension(ndims), intent(in) :: periods
      logical :: reorder = .TRUE.

      integer(kind=4), intent(out) :: comm_cart, ierr

      call MPI_Cart_create(INT(MPI_COMM_WORLD,kind=4), INT(ndims,kind=4), INT(dims,kind=4),&
         periods, reorder, comm_cart, ierr)

   end subroutine create_cart_communicator_mpi_wrapper

   !> Get the Cartesian coordinates of the current process
   subroutine get_cart_coords_mpi_wrapper(comm, rank, ndims, coords, ierr)
      integer, intent(in) :: comm, rank, ndims

      integer(kind=4), intent(out) :: coords(ndims), ierr

      call MPI_Cart_coords(INT(comm,kind=4), INT(rank,kind=4), INT(ndims,kind=4), coords, ierr)
   end subroutine get_cart_coords_mpi_wrapper

   !> Routine to get the rank of a process given its coordinates
   subroutine mpi_cart_rank_mpi_wrapper(comm, ndims, indices, rank_of_coords, error, ierr)
      integer, intent(in) :: comm, ndims
      integer, dimension(ndims), intent(in) :: indices

      integer(kind=4), intent(out) :: rank_of_coords, error, ierr

      error = 0
      call MPI_Cart_rank(INT(comm,kind=4), INT(indices,kind=4), rank_of_coords, ierr)
      if(ierr /= MPI_SUCCESS) then
         error = 1
      endif
   end subroutine mpi_cart_rank_mpi_wrapper

   !> Routine to change the MPI_COMM error handler so we can find neighbors without crashing. We restore the original using original_MPI_COMM_errhandler() when done
   subroutine change_MPI_COMM_errhandler_mpi_wrapper(comm, ierr)
      integer, intent(in) :: comm

      integer(kind=4), intent(out) :: ierr

      ! Check if the original error handler has already been saved
      if (original_errhandler == MPI_ERRHANDLER_NULL) then
         ! Get the current error handler for comm
         call MPI_Comm_get_errhandler(INT(comm,kind=4), original_errhandler, ierr)
         if (ierr /= MPI_SUCCESS) then
            print *, "Error getting original MPI error handler."
            return
         endif
      endif

      ! Set the error handler for comm to return errors
      call MPI_Comm_set_errhandler(INT(comm,kind=4), MPI_ERRORS_RETURN, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, "Error setting MPI error handler to MPI_ERRORS_RETURN."
      endif
   end subroutine change_MPI_COMM_errhandler_mpi_wrapper

   !> Restores the MPI_COMM error handler back to the original after finding neighbors
   subroutine original_MPI_COMM_errhandler_mpi_wrapper(comm, ierr)
      integer, intent(in) :: comm

      integer(kind=4), intent(out) :: ierr

      ! Check if the original error handler was saved
      if (original_errhandler /= MPI_ERRHANDLER_NULL) then
         ! Restore the original error handler for parameters%cart_comm
         call MPI_Comm_set_errhandler(INT(comm,kind=4), original_errhandler, ierr)
         if (ierr /= MPI_SUCCESS) then
            print *, "Error restoring original MPI error handler."
         endif
         ! Reset the original_errhandler variable
         original_errhandler = MPI_ERRHANDLER_NULL
      else
         print *, "Original MPI error handler not saved."
      endif
   end subroutine original_MPI_COMM_errhandler_mpi_wrapper

end module mpi_wrapper_module
