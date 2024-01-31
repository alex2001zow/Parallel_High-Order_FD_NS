module rank_parameters_module
   use mpi
   use utility_functions_module, only : IDX_XD
   use constants_module, only: neighbor_current_rank, neighbor_non_existant_rank
   implicit none

   private
   !character(len=MPI_MAX_ERROR_STRING) :: error_string
   !integer :: ierr_string, error_len

   integer :: original_errhandler = MPI_ERRHANDLER_NULL
   integer :: ierr = 0

   !> Structure to hold the parameters for each rank.
   type rank_struct
      integer :: ndims, total_num_elements
      integer :: rank, world_size, cart_comm
      integer :: total_number_of_neighbors, effective_number_of_neighbors
      integer, allocatable :: grid_size(:), processor_dim(:), coords(:), neighbors(:), local_size(:)
      logical, allocatable :: periods(:)
   end type rank_struct

   public :: rank_struct, setup_rank_parameters, deallocate_rank_parameters, print_rank_parameters

contains

   !> Allocate and setup each ranks parameters
   subroutine setup_rank_parameters(parameters, rank_val, world_size_val)
      character(255) :: dim_input_arg, temp_arg
      integer, intent(in) :: rank_val, world_size_val
      type(rank_struct), intent(out) :: parameters

      integer :: ii

      call get_command_argument(1, dim_input_arg)
      read(dim_input_arg,*) parameters%ndims

      parameters%rank = rank_val
      parameters%world_size = world_size_val

      parameters%total_number_of_neighbors = 3**parameters%ndims

      call allocate_rank_parameters(parameters)

      parameters%total_num_elements = 1
      do ii = 1, parameters%ndims
         call get_command_argument(ii+1, temp_arg)
         read(temp_arg,*) parameters%grid_size(ii)

         call get_command_argument(ii+1+parameters%ndims, temp_arg)
         read(temp_arg,*) parameters%processor_dim(ii)

         ! We have make sure that this is properly divisible
         if(mod(parameters%grid_size(ii), parameters%processor_dim(ii)) /= 0) then
            print *, "Grid size is not divisible by the number of processors in dimension ", ii
            stop
         end if
         parameters%local_size(ii) = parameters%grid_size(ii) / parameters%processor_dim(ii)

         ! We do not have any period in any of the dimensions in the cartesian processor system
         parameters%periods(ii) = .FALSE.

         ! The total number of elements we have to allocate per node
         parameters%total_num_elements = parameters%total_num_elements * parameters%local_size(ii)

      end do

      ! Create a Cartesian communicator
      call MPI_Cart_create(MPI_COMM_WORLD, parameters%ndims, parameters%processor_dim, parameters%periods, .TRUE., &
         parameters%cart_comm, ierr)

      ! Get the Cartesian coordinates of the current process
      call MPI_Cart_coords(parameters%cart_comm, parameters%rank, parameters%ndims, parameters%coords, ierr)

      if(parameters%ndims == 2) then
         call determine_moore_neighbors_2d(parameters)
      else if(parameters%ndims == 3) then
         call determine_moore_neighbors_3d(parameters)
      else
         print *, "Only 2D and 3D is currently supported"
         stop
      endif

   end subroutine setup_rank_parameters

   !> Allocate rank parameters. IMPORTANT neighbors are not the correct size yet!
   subroutine allocate_rank_parameters(parameters)
      type(rank_struct), intent(inout) :: parameters

      allocate(parameters%grid_size(parameters%ndims))
      allocate(parameters%processor_dim(parameters%ndims))
      allocate(parameters%periods(parameters%ndims))
      allocate(parameters%coords(parameters%ndims))
      allocate(parameters%neighbors(parameters%total_number_of_neighbors))
      allocate(parameters%local_size(parameters%ndims))
   end subroutine allocate_rank_parameters

   !> Deallocate rank parameters
   subroutine deallocate_rank_parameters(parameters)
      type(rank_struct), intent(inout) :: parameters

      ! Deallocate the pointers
      deallocate(parameters%grid_size)
      deallocate(parameters%processor_dim)
      deallocate(parameters%periods)
      deallocate(parameters%coords)
      deallocate(parameters%neighbors)
      deallocate(parameters%local_size)

      call MPI_Comm_free(parameters%cart_comm, ierr)

   end subroutine deallocate_rank_parameters

   !> Print rank parameters
   subroutine print_rank_parameters(parameters)
      type(rank_struct), intent(in) :: parameters
      print *, "Dimensions:", parameters%ndims
      print *, "Rank:", parameters%rank
      print *, "World Size:", parameters%world_size
      print *, "Grid Size:", parameters%grid_size
      print *, "Num Processors per dim:", parameters%processor_dim
      print *, "Periodic:", parameters%periods
      print *, "Cartesian coordinates:", parameters%coords
      print *, "Neighbors:", parameters%neighbors
      print *, "Local matrix size:", parameters%local_size
      print *, "Total number of elements:", parameters%total_num_elements
   end subroutine print_rank_parameters

   !> Find the neighbors of each rank in a 2D grid. Would like to make it N-D but that is a bit more complicated.
   subroutine determine_moore_neighbors_2d(parameters)
      type(rank_struct), intent(inout) :: parameters

      integer :: ii, jj, global_index, rank_of_coords
      integer :: indices(parameters%ndims)

      ! Prevent MPI from aborting when we try to find neighbors that do not exist
      call change_MPI_COMM_errhandler(parameters)

      ! Determine the ranks of neighboring processes including corners along each dimension.
      global_index = 1
      do ii = -1,1,1
         do jj = -1,1,1
            indices = parameters%coords + [ii,jj]

            call MPI_Cart_rank(parameters%cart_comm, indices, rank_of_coords, ierr)

            if(ierr /= MPI_SUCCESS) then
               parameters%neighbors(global_index) = neighbor_non_existant_rank ! This is a non-existent rank
            else if(ii == 0 .AND. jj == 0) then
               parameters%neighbors(global_index) = neighbor_current_rank ! This is the current rank
            else
               parameters%neighbors(global_index) = rank_of_coords
            endif
            global_index = global_index + 1
         end do
      end do

      ! Restore the original error handler to abort on errors
      call original_MPI_COMM_errhandler(parameters)

   end subroutine determine_moore_neighbors_2d

   !> Find the neighbors of each rank in a 3D grid. Would like to make it N-D but that is a bit more complicated.
   subroutine determine_moore_neighbors_3d(parameters)
      type(rank_struct), intent(inout) :: parameters

      integer :: ii, jj, kk, global_index, rank_of_coords
      integer :: indices(parameters%ndims)

      ! Prevent MPI from aborting when we try to find neighbors that do not exist
      call change_MPI_COMM_errhandler(parameters)

      ! Determine the ranks of neighboring processes including corners along each dimension.
      global_index = 1
      do ii = -1,1,1
         do jj = -1,1,1
            do kk = -1,1,1
               indices = parameters%coords + [ii,jj,kk]

               call MPI_Cart_rank(parameters%cart_comm, indices, rank_of_coords, ierr)

               if(ierr /= MPI_SUCCESS) then
                  parameters%neighbors(global_index) = neighbor_non_existant_rank ! This is a non-existent rank
               else if(ii == 0 .AND. jj == 0 .AND. kk == 0) then
                  parameters%neighbors(global_index) = neighbor_current_rank ! This is the current rank
               else
                  parameters%neighbors(global_index) = rank_of_coords
               endif
               global_index = global_index + 1
            end do
         end do
      end do

      ! Restore the original error handler to abort on errors
      call original_MPI_COMM_errhandler(parameters)

   end subroutine determine_moore_neighbors_3d

   !> Routine to find the elements to be sent and recieved from each neighbor
   !> todo

   !> Routine to send and reviece elements from each neighbor
   !> todo

   !> Routine to change the MPI_COMM error handler so we can find neighbors without crashing. We restore the original using original_MPI_COMM_errhandler() when done
   subroutine change_MPI_COMM_errhandler(parameters)
      type(rank_struct), intent(inout) :: parameters

      ! Check if the original error handler has already been saved
      if (original_errhandler == MPI_ERRHANDLER_NULL) then
         ! Get the current error handler for parameters%cart_comm
         call MPI_Comm_get_errhandler(parameters%cart_comm, original_errhandler, ierr)
         if (ierr /= MPI_SUCCESS) then
            print *, "Error getting original MPI error handler."
            return
         endif
      endif

      ! Set the error handler for parameters%cart_comm to return errors
      call MPI_Comm_set_errhandler(parameters%cart_comm, MPI_ERRORS_RETURN, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, "Error setting MPI error handler to MPI_ERRORS_RETURN."
      endif
   end subroutine change_MPI_COMM_errhandler

   !> Restores the MPI_COMM error handler back to the original after finding neighbors
   subroutine original_MPI_COMM_errhandler(parameters)
      type(rank_struct), intent(inout) :: parameters

      ! Check if the original error handler was saved
      if (original_errhandler /= MPI_ERRHANDLER_NULL) then
         ! Restore the original error handler for parameters%cart_comm
         call MPI_Comm_set_errhandler(parameters%cart_comm, original_errhandler, ierr)
         if (ierr /= MPI_SUCCESS) then
            print *, "Error restoring original MPI error handler."
         endif
         ! Reset the original_errhandler variable
         original_errhandler = MPI_ERRHANDLER_NULL
      else
         print *, "Original MPI error handler not saved."
      endif
   end subroutine original_MPI_COMM_errhandler
end module rank_parameters_module
