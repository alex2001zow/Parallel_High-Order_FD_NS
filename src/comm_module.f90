module comm_module
   use mpi, only : MPI_REQUEST_NULL
   use constants_module, only: neighbor_current_rank, neighbor_non_existant_rank
   use utility_functions_module, only: IDX_XD, IDX_XD_INV, print_integer_1D_array, print_integer_2D_array, &
      sleeper_function, reshape_integer_1D_to_2D, reshape_integer_1D_to_3D
   use mpi_wrapper_module, only: create_cart_communicator_mpi_wrapper, &
      get_cart_coords_mpi_wrapper,change_MPI_COMM_errhandler_mpi_wrapper, &
      original_MPI_COMM_errhandler_mpi_wrapper, cart_rank_mpi_wrapper, &
      free_communicator_mpi_wrapper
   implicit none

   private

   type comm_type
      integer :: ndims, rank, world_size, comm
      logical :: reorder
      integer, dimension(:), allocatable :: processor_dim, coords
      logical, dimension(:), allocatable :: periods
   end type comm_type

   public :: comm_type, create_cart_comm_type, deallocate_cart_comm_type, print_cart_comm_type, &
      get_neighbor_indices
   public :: get_cart_neighbors, begin_end_neighbor_indices

contains

   !> Allocate a comm_type object
   subroutine create_cart_comm_type(ndims, processor_dim, periods, reorder, rank, world_size, comm_output)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(:), intent(in) :: processor_dim
      logical, dimension(:), intent(in) :: periods
      logical, intent(in) :: reorder
      type(comm_type), intent(out) :: comm_output

      comm_output%ndims = ndims
      comm_output%rank = rank
      comm_output%world_size = world_size

      call allocate_cart_comm_type(ndims, comm_output)

      comm_output%processor_dim = processor_dim
      comm_output%periods = periods
      comm_output%reorder = reorder

      call create_cart_communicator(comm_output%ndims, comm_output%processor_dim, comm_output%periods, &
         comm_output%reorder, comm_output%comm)

      call get_cart_coords(comm_output%comm, rank, ndims, comm_output%coords)

   end subroutine create_cart_comm_type

   !> Create a cartesian communicator in OpenMPI
   subroutine create_cart_communicator(ndims, processor_dim, periods, reorder, cart_comm)
      integer, intent(in) :: ndims
      logical, intent(in) :: reorder
      integer, dimension(:), intent(in) :: processor_dim
      logical, dimension(:), intent(in) :: periods

      integer, intent(out) :: cart_comm

      call create_cart_communicator_mpi_wrapper(ndims, processor_dim, periods, reorder, cart_comm)

   end subroutine create_cart_communicator

   !> Get the cartesian coordinates of a rank in a communicator
   subroutine get_cart_coords(cart_comm, rank, ndims, coords)
      integer, intent(in) :: cart_comm, rank, ndims
      integer, dimension(:), intent(out) :: coords

      call get_cart_coords_mpi_wrapper(cart_comm, rank, ndims, coords)

   end subroutine get_cart_coords

   !> Find the Moore neighborhood of a rank in a cartesian communicator. Gives out the neighbors, ghost points, and offsets. Ghost points and offsets are always positive.
   !! If they 1 then we have a boundary between the blocks or the true boundary.
   subroutine get_cart_neighbors(ndims, coords, cart_comm, neighbors, &
      begin_ghost_offset, end_ghost_offset, begin_stencil_offset, end_stencil_offset)
      integer, intent(in) :: ndims, cart_comm
      integer, dimension(:), intent(in) :: coords
      integer, dimension(:), intent(out) :: neighbors
      integer, dimension(:), intent(out) :: begin_ghost_offset, end_ghost_offset, begin_stencil_offset, end_stencil_offset

      integer :: global_index, num_nonzero
      integer, dimension(ndims) :: neighbor_array, indices, current_index
      integer :: rank_of_coords, error

      ! Prevent MPI from aborting when we try to find neighbors that do not exist
      call change_MPI_COMM_errhandler_mpi_wrapper(cart_comm)

      ! Determine the ghost points for the current rank
      begin_ghost_offset = 0
      end_ghost_offset = 0

      ! Determine the ranks of neighboring processes including corners along each dimension.
      begin_stencil_offset = 0
      end_stencil_offset = 0

      neighbor_array = 3 ! This is the number of neighbors in each dimension

      ! Run through all the neighbors
      do global_index = 1, 3**ndims
         call IDX_XD_INV(ndims, neighbor_array, global_index, current_index)
         current_index = current_index - 2 ! To go -1, 0, 1

         indices = coords + current_index

         call cart_rank_mpi_wrapper(cart_comm, ndims, indices, rank_of_coords, error)

         num_nonzero = count(current_index /= 0)  ! Count the number of non-zero entries in current_index

         ! error = 1 means no neighbor exists. We have true boundary.
         if(error == 1) then
            neighbors(global_index) = neighbor_non_existant_rank ! This is a non-existent rank. Should be replaced with a MPI_PROC_NULL. Do it in the constants module

            if(ndims == 1 .or. num_nonzero /= ndims) then
               begin_ghost_offset = min(begin_ghost_offset, current_index) ! Multiply with ghost size. Done elsewhere. Seems to be correct
               end_ghost_offset = max(end_ghost_offset, current_index) ! Multiply with ghost size. Done elsewhere. Seems to be correct
            end if
         else if(all(current_index == 0)) then
            neighbors(global_index) = neighbor_current_rank ! This is the current rank
         else
            neighbors(global_index) = rank_of_coords ! This is a neighbor rank
            begin_stencil_offset = min(begin_stencil_offset, current_index) ! Multiply with stencil size. Done elsewhere. This is correct
            end_stencil_offset = max(end_stencil_offset, current_index) ! Multiply with stencil size. Done elsewhere. This is correct
         endif

      end do

      ! Take the absolute value of the points
      begin_ghost_offset = abs(begin_ghost_offset)
      end_ghost_offset = abs(end_ghost_offset)
      begin_stencil_offset = abs(begin_stencil_offset)
      end_stencil_offset = abs(end_stencil_offset)

      ! Restore the original error handler to abort on errors
      call original_MPI_COMM_errhandler_mpi_wrapper(cart_comm)

   end subroutine get_cart_neighbors

   !> Calculate the indices of the elements that need to be sent to the neighbors from the block. Can take multiple dimensions. send_or_recv = 0 for send, 1 for recv.
   pure subroutine begin_end_neighbor_indices(ndims, extended_block_begin_c, extended_block_end_c, block_begin_c, block_end_c, &
      stencil_size, indices, send_or_recv, begin, end, output_dims)
      integer, intent(in) :: ndims, send_or_recv
      integer, dimension(:), intent(in) :: extended_block_begin_c, extended_block_end_c
      integer, dimension(:), intent(in) :: block_begin_c, block_end_c
      integer, dimension(:), intent(in) :: stencil_size, indices
      integer, dimension(:), intent(out):: begin, end, output_dims

      integer :: current_dim

      ! For the send indices
      do current_dim = 1, ndims
         select case(send_or_recv)
          case(0)
            select case(indices(current_dim))
             case(-1)
               begin(current_dim) = block_begin_c(current_dim)
               end(current_dim) = block_begin_c(current_dim) + stencil_size(current_dim)
             case(0)
               begin(current_dim) = block_begin_c(current_dim)
               end(current_dim) = block_end_c(current_dim)
             case(1)
               begin(current_dim) = block_end_c(current_dim) - stencil_size(current_dim)
               end(current_dim) = block_end_c(current_dim)
            end select

            ! For the recv indices
          case(1)
            select case(indices(current_dim))
             case(-1)
               begin(current_dim) = extended_block_begin_c(current_dim)
               end(current_dim) = extended_block_begin_c(current_dim) + block_begin_c(current_dim)
             case(0)
               begin(current_dim) = block_begin_c(current_dim)
               end(current_dim) = block_end_c(current_dim)
             case(1)
               begin(current_dim) = block_end_c(current_dim)
               end(current_dim) = extended_block_end_c(current_dim)
            end select
         end select
      end do

      output_dims = end - begin + 1

   end subroutine begin_end_neighbor_indices

   !> Subroutine to get the neighbor indices out
   pure subroutine get_neighbor_indices(ndims, neighbor_index, indices, begin, end)
      integer, intent(in) :: ndims, neighbor_index
      integer, dimension(:), intent(in) :: indices
      integer, dimension(:), intent(out) :: begin, end

      integer :: start_index

      start_index = (neighbor_index-1) * 2*ndims + 1

      begin = indices(start_index: start_index + ndims - 1)
      end = indices(start_index + ndims: start_index + 2*ndims - 1)

   end subroutine get_neighbor_indices

   !> Allocate a comm_type object
   pure subroutine allocate_cart_comm_type(ndims, comm)
      integer, intent(in) :: ndims
      type(comm_type), intent(inout) :: comm

      allocate(comm%processor_dim(ndims))
      allocate(comm%periods(ndims))
      allocate(comm%coords(ndims))

   end subroutine allocate_cart_comm_type

   !> Deallocate a comm_type object
   subroutine deallocate_cart_comm_type(comm)
      type(comm_type), intent(inout) :: comm

      call free_communicator_mpi_wrapper(comm%comm)

      if(allocated(comm%processor_dim)) deallocate(comm%processor_dim)
      if(allocated(comm%periods)) deallocate(comm%periods)
      if(allocated(comm%coords)) deallocate(comm%coords)

   end subroutine deallocate_cart_comm_type

   !> Print the contents of a comm_type object
   subroutine print_cart_comm_type(comm, iounit)
      type(comm_type), intent(in) :: comm
      integer, intent(in) :: iounit

      integer, dimension(comm%ndims) :: neighbor_array
      integer :: current_rank, start_index, end_index
      integer, dimension(comm%ndims) :: coords
      integer, dimension(comm%world_size), target :: ranks
      integer, dimension(comm%world_size * (comm%ndims+1)), target :: rank_and_coords
      integer, contiguous, dimension(:), pointer :: ranks_ptr
      integer, contiguous, dimension(:,:), pointer :: ranks_2D_ptr

      neighbor_array = 3

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "comm_type: "
      write(iounit, *)

      write(iounit, *) "ndims: ", comm%ndims
      write(iounit, *) "rank: ", comm%rank
      write(iounit, *) "world_size: ", comm%world_size
      write(iounit, *) "comm: ", comm%comm
      write(iounit, *) "reorder: ", comm%reorder

      write(iounit, *) "processor_dim: ", comm%processor_dim

      write(iounit, *) "coords: ", comm%coords
      write(iounit, *) "periods: ", comm%periods

      write(iounit, *) "Cartesian grid:"

      ! Print the cartesian grid of the ranks.

      ! Iterate over all ranks in the cartesian communicator
      do current_rank = 1, comm%world_size
         ! Get the coordinates for the current rank
         call get_cart_coords_mpi_wrapper(comm%comm, current_rank-1, comm%ndims, coords)
         ranks(current_rank) = current_rank-1

         start_index = (current_rank-1) * (comm%ndims+1) + 1
         end_index = current_rank * (comm%ndims+1)
         rank_and_coords(start_index:end_index) = [current_rank-1, coords]

      end do

      ranks_ptr => ranks

      select case(comm%ndims)
       case(1)
         call print_integer_1D_array(ranks, iounit)
       case(2)
         call reshape_integer_1D_to_2D(comm%processor_dim, ranks, ranks_2D_ptr)
         call print_integer_2D_array(ranks_2D_ptr, iounit)
       case default
         write(iounit, *) "Printing for more than 2 dimensions is not supported."
      end select

   end subroutine print_cart_comm_type

end module comm_module
