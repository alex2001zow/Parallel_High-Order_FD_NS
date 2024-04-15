module comm_module
   use mpi, only : MPI_REQUEST_NULL
   use constants_module, only: neighbor_current_rank, neighbor_non_existant_rank
   use utility_functions_module, only: IDX_XD, IDX_XD_INV, print_real_array, print_integer_array, sleeper_function
   use mpi_wrapper_module, only: create_cart_communicator_mpi_wrapper, &
      get_cart_coords_mpi_wrapper,change_MPI_COMM_errhandler_mpi_wrapper, &
      original_MPI_COMM_errhandler_mpi_wrapper, cart_rank_mpi_wrapper, &
      free_communicator_mpi_wrapper
   implicit none

   private

   type comm_type
      integer :: comm, num_neighbors
      integer, dimension(:), allocatable :: coords, neighbors, offset_begin, offset_end
      integer, dimension(:), allocatable :: begin, end, size
      integer, dimension(:), allocatable :: neighbor_send_indices, neighbor_recv_indices, neighbor_elements_size
      integer, dimension(:), allocatable :: neighbor_send_request, neighbor_recv_request
   end type comm_type

   public :: comm_type, create_cart_comm_type, deallocate_cart_comm_type, print_cart_comm_type, print_cartesian_grid, &
      get_neighbor_indices

contains

   !> Allocate a comm_type object
   subroutine create_cart_comm_type(ndims, dimensions, processor_dim, rank, comm_output)
      integer, intent(in) :: ndims, rank
      integer, dimension(:), intent(in) :: dimensions, processor_dim
      type(comm_type), intent(inout) :: comm_output

      comm_output%num_neighbors = 3**ndims

      call allocate_cart_comm_type(ndims, comm_output)

      call create_cart_communicator(ndims, processor_dim, comm_output%comm)

      call get_cart_coords(comm_output%comm, rank, ndims, comm_output%coords)

      call get_cart_neighbors(ndims, comm_output%coords, comm_output%comm, comm_output%neighbors, &
         comm_output%offset_begin, comm_output%offset_end)

      comm_output%size = (dimensions/processor_dim)

      comm_output%begin = (comm_output%coords * comm_output%size) + 1 + comm_output%offset_begin
      comm_output%end = (comm_output%begin - comm_output%offset_begin) + comm_output%size - 1 + comm_output%offset_end

      comm_output%size = comm_output%end - comm_output%begin + 1

      call calculate_neighbor_indices(ndims, comm_output%size, comm_output%neighbor_send_indices, comm_output%neighbor_recv_indices)

      call calculate_number_of_elements_to_send(ndims, comm_output)

      comm_output%neighbor_send_request(:) = MPI_REQUEST_NULL
      comm_output%neighbor_recv_request(:) = MPI_REQUEST_NULL

   end subroutine create_cart_comm_type

   !> Create a cartesian communicator in OpenMPI
   subroutine create_cart_communicator(ndims, processor_dim, cart_comm)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: processor_dim

      integer, intent(out) :: cart_comm

      call create_cart_communicator_mpi_wrapper(ndims, processor_dim, cart_comm)

   end subroutine create_cart_communicator

   !> Get the cartesian coordinates of a rank in a communicator
   subroutine get_cart_coords(cart_comm, rank, ndims, coords)
      integer, intent(in) :: cart_comm, rank, ndims
      integer, dimension(:), intent(inout) :: coords

      call get_cart_coords_mpi_wrapper(cart_comm, rank, ndims, coords)

   end subroutine get_cart_coords

   !> Find the Moore neighborhood of a rank in a cartesian communicator. We should check the order column vs row. Not sure it is correct but it seems to work.
   subroutine get_cart_neighbors(ndims, coords, cart_comm, neighbors, offset_begin, offset_end)
      integer, intent(in) :: ndims, cart_comm
      integer, dimension(:), intent(in) :: coords
      integer, dimension(:), intent(inout) :: neighbors
      integer, dimension(:), intent(inout) :: offset_begin, offset_end

      integer :: global_index
      integer, dimension(ndims) :: neighbor_array, indices, current_index
      integer :: rank_of_coords, error

      ! Prevent MPI from aborting when we try to find neighbors that do not exist
      call change_MPI_COMM_errhandler_mpi_wrapper(cart_comm)

      ! Determine the ranks of neighboring processes including corners along each dimension.
      offset_begin(:) = 0
      offset_end(:) = 0

      neighbor_array(:) = 3 ! This is the number of neighbors in each dimension

      do global_index = 1, 3**ndims
         call IDX_XD_INV(ndims, neighbor_array, global_index, current_index)
         current_index = current_index - 2 ! To go -1, 0, 1

         indices = coords + current_index

         call cart_rank_mpi_wrapper(cart_comm, ndims, indices, rank_of_coords, error)

         if(error == 1) then
            neighbors(global_index) = neighbor_non_existant_rank ! This is a non-existent rank
         else if(all(current_index == 0)) then
            neighbors(global_index) = neighbor_current_rank ! This is the current rank
         else
            neighbors(global_index) = rank_of_coords
            offset_begin(:) = min(offset_begin(:), current_index)
            offset_end(:) = max(offset_end(:), current_index)
         endif
      end do

      ! Restore the original error handler to abort on errors
      call original_MPI_COMM_errhandler_mpi_wrapper(cart_comm)

   end subroutine get_cart_neighbors

   !> Calculate the indices of the elements that need to be sent to the neighbors
   pure subroutine calculate_neighbor_indices(ndims, dims, send_array, recv_array)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: dims
      integer, dimension(:), intent(inout) :: send_array, recv_array

      integer :: global_index, start_index, end_index
      integer, dimension(2*ndims) :: send_indices, recv_indices ! begin and end for each dimension
      integer, dimension(ndims) :: neighbor_array, current_index

      neighbor_array(:) = 3

      send_array(:) = -10
      recv_array(:) = -10

      global_index = 1

      do global_index = 1, 3**ndims
         call IDX_XD_INV(ndims, neighbor_array, global_index, current_index)
         current_index = current_index - 2 ! To go -1, 0, 1

         start_index = (global_index-1) * 2*ndims + 1
         end_index = global_index * 2*ndims

         ! Send indices
         call begin_end_neighbor_indices(ndims, dims, current_index, 0, send_indices(1:ndims), send_indices(ndims + 1:2*ndims))
         ! Recv indices
         call begin_end_neighbor_indices(ndims, dims, current_index, 1, recv_indices(1:ndims), recv_indices(ndims + 1:2*ndims))

         send_array(start_index:end_index) = send_indices(:)
         recv_array(start_index:end_index) = recv_indices(:)

      end do

   end subroutine calculate_neighbor_indices

   !> Calculate the indices of the elements that need to be sent to the neighbors from the block. Can take multiple dimensions. send_or_recv = 0 for send, 1 for recv
   pure subroutine begin_end_neighbor_indices(ndims, dims, indices, send_or_recv, begin, end)
      integer, intent(in) :: ndims, send_or_recv
      integer, dimension(:), intent(in) :: dims, indices
      integer, dimension(:), intent(out):: begin, end

      integer :: current_dim

      ! For the send indices
      if(send_or_recv == 0) then
         do current_dim = 1, ndims
            if(indices(current_dim) == -1) then
               begin(current_dim) = 2
               end(current_dim) = 2
            else if(indices(current_dim) == 0) then
               begin(current_dim)  = 2
               end(current_dim)  = dims(current_dim)  - 1
            else if(indices(current_dim) == 1) then
               begin(current_dim)  = dims(current_dim)  - 1
               end(current_dim)  = dims(current_dim)  - 1
            endif
         end do
      end if

      ! For the recv indices
      if(send_or_recv == 1) then
         do current_dim = 1, ndims
            if(indices(current_dim) == -1) then
               begin(current_dim) = 1
               end(current_dim) = 1
            else if(indices(current_dim) == 0) then
               begin(current_dim)  = 2
               end(current_dim)  = dims(current_dim)  - 1
            else if(indices(current_dim) == 1) then
               begin(current_dim)  = dims(current_dim)
               end(current_dim)  = dims(current_dim)
            endif
         end do
      end if


   end subroutine begin_end_neighbor_indices

   !> Subroutine to get the neighbor indices out
   pure subroutine get_neighbor_indices(ndims, neighbor_index, indices, begin, end)
      integer, intent(in) :: ndims, neighbor_index
      integer, dimension(:), intent(in) :: indices
      integer, dimension(:), intent(out) :: begin, end

      integer :: start_index

      start_index = (neighbor_index-1) * 2*ndims + 1

      begin(:) = indices(start_index: start_index + ndims - 1)
      end(:) = indices(start_index + ndims: start_index + 2*ndims - 1)

   end subroutine get_neighbor_indices

   !> Subroutine to calculate the number of elements that need to be sent to the neighbors
   pure subroutine calculate_number_of_elements_to_send(ndims, comm_inout)
      integer, intent(in) :: ndims
      type(comm_type), intent(inout) :: comm_inout

      integer :: neighbor_index
      integer, dimension(ndims) :: begin, end

      do neighbor_index = 1, comm_inout%num_neighbors
         call get_neighbor_indices(ndims, neighbor_index, comm_inout%neighbor_send_indices, begin, end)
         comm_inout%neighbor_elements_size(neighbor_index) = product(end - begin + 1)
      end do

   end subroutine calculate_number_of_elements_to_send

   !> Allocate a comm_type object
   pure subroutine allocate_cart_comm_type(ndims, comm)
      integer, intent(in) :: ndims
      type(comm_type), intent(inout) :: comm

      allocate(comm%coords(ndims))
      allocate(comm%neighbors(comm%num_neighbors))
      allocate(comm%offset_begin(ndims))
      allocate(comm%offset_end(ndims))
      allocate(comm%begin(ndims))
      allocate(comm%end(ndims))
      allocate(comm%size(ndims))
      allocate(comm%neighbor_send_indices(comm%num_neighbors*2*ndims))
      allocate(comm%neighbor_recv_indices(comm%num_neighbors*2*ndims))
      allocate(comm%neighbor_elements_size(comm%num_neighbors))
      allocate(comm%neighbor_send_request(comm%num_neighbors))
      allocate(comm%neighbor_recv_request(comm%num_neighbors))

   end subroutine allocate_cart_comm_type

   !> Deallocate a comm_type object
   subroutine deallocate_cart_comm_type(comm)
      type(comm_type), intent(inout) :: comm

      call free_communicator_mpi_wrapper(comm%comm)

      if(allocated(comm%coords)) deallocate(comm%coords)
      if(allocated(comm%neighbors)) deallocate(comm%neighbors)
      if(allocated(comm%offset_begin)) deallocate(comm%offset_begin)
      if(allocated(comm%offset_end)) deallocate(comm%offset_end)
      if(allocated(comm%begin)) deallocate(comm%begin)
      if(allocated(comm%end)) deallocate(comm%end)
      if(allocated(comm%size)) deallocate(comm%size)
      if(allocated(comm%neighbor_send_indices)) deallocate(comm%neighbor_send_indices)
      if(allocated(comm%neighbor_recv_indices)) deallocate(comm%neighbor_recv_indices)
      if(allocated(comm%neighbor_elements_size)) deallocate(comm%neighbor_elements_size)
      if(allocated(comm%neighbor_send_request)) deallocate(comm%neighbor_send_request)
      if(allocated(comm%neighbor_recv_request)) deallocate(comm%neighbor_recv_request)

   end subroutine deallocate_cart_comm_type

   !> Print the contents of a comm_type object
   subroutine print_cart_comm_type(ndims, comm, iounit)
      integer, intent(in) :: ndims
      type(comm_type), intent(in) :: comm
      integer, intent(in) :: iounit

      integer, dimension(ndims) :: neighbor_array

      neighbor_array(:) = 3

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "comm_type: "
      write(iounit, *)

      write(iounit, *) "comm: ", comm%comm
      write(iounit, *) "coords: ", comm%coords

      write(iounit, *) "begin: ", comm%begin
      write(iounit, *) "end: ", comm%end
      write(iounit, *) "offset_begin: ", comm%offset_begin
      write(iounit, *) "offset_end: ", comm%offset_end
      write(iounit, *) "size: ", comm%size

      call print_integer_array(ndims, neighbor_array, comm%neighbor_elements_size, 1, "Neighbor element size: ", iounit)

      call print_integer_array(ndims, neighbor_array, comm%neighbor_send_indices, 2*ndims, "Neighbor send indices: ", iounit)

      call print_integer_array(ndims, neighbor_array, comm%neighbor_recv_indices, 2*ndims, "Neighbor recv indices: ", iounit)

      call print_integer_array(ndims, neighbor_array, comm%neighbors, 1, "Neighbor ranks: ", iounit)

      call print_integer_array(ndims, neighbor_array, comm%neighbor_send_request, 1, "Neighbor send request: ", iounit)

      call print_integer_array(ndims, neighbor_array, comm%neighbor_recv_request, 1, "Neighbor recv request: ", iounit)

   end subroutine print_cart_comm_type

   !> A routine to print the cartesian grid of the ranks.
   subroutine print_cartesian_grid(ndims, cart_comm, world_size, processor_dim, filename)
      integer, intent(in) :: ndims, cart_comm, world_size
      integer, dimension(:), intent(in) :: processor_dim
      character(len=*), intent(in) :: filename

      integer :: current_rank, start_index, end_index, iounit, ios
      integer, dimension(ndims) :: coords
      integer, dimension(world_size * (ndims+1)) :: rank_and_coords
      character(255) :: file_with_grid

      ! Construct the filename by appending "_cart_grid" to the original filename
      write(file_with_grid, '(A, "cart_grid.txt")') trim(filename)

      ! Open the file for writing, replacing it if it already exists
      open(newunit=iounit, file=trim(file_with_grid), status='replace', action='write', iostat=ios)

      ! Check for errors in opening the file
      if (ios /= 0) then
         print *, "Error opening file: ", trim(file_with_grid)
         return
      endif

      ! Write the header to the file
      write(iounit, *) "Rank cartesian processor grid with dimension:", processor_dim

      ! Iterate over all ranks in the cartesian communicator
      do current_rank = 1, world_size
         ! Get the coordinates for the current rank
         call get_cart_coords_mpi_wrapper(cart_comm, current_rank-1, ndims, coords)

         start_index = (current_rank-1) * (ndims+1) + 1
         end_index = current_rank * (ndims+1)
         rank_and_coords(start_index:end_index) = [current_rank-1, coords]

      end do

      call print_integer_array(ndims, processor_dim, rank_and_coords, ndims+1, "Rank and coords: ", iounit)

      ! Close the file
      close(iounit)
   end subroutine print_cartesian_grid

end module comm_module
