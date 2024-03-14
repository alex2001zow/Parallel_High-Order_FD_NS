module comm_module
   use mpi, only : MPI_REQUEST_NULL
   use constants_module, only: neighbor_current_rank, neighbor_non_existant_rank
   use utility_functions_module, only: IDX_XD
   use mpi_wrapper_module, only: create_cart_communicator_mpi_wrapper, &
      get_cart_coords_mpi_wrapper,change_MPI_COMM_errhandler_mpi_wrapper, &
      original_MPI_COMM_errhandler_mpi_wrapper, cart_rank_mpi_wrapper, &
      free_communicator_mpi_wrapper
   implicit none

   private

   type comm_type
      integer :: comm, num_neighbors
      integer, allocatable :: coords(:), neighbors(:), offset_begin(:), offset_end(:)
      integer, allocatable :: begin(:), end(:), size(:)
      integer, allocatable :: neighbor_send_indices(:), neighbor_recv_indices(:), neighbor_elements_size(:)
      integer, allocatable :: neighbor_send_request(:), neighbor_recv_request(:)
   end type comm_type

   public :: comm_type, create_cart_comm_type, deallocate_cart_comm_type, print_cart_comm_type, print_cartesian_grid, &
      get_neighbor_indices

contains

   !> Allocate a comm_type object
   subroutine create_cart_comm_type(ndims, dimensions, processor_dim, rank, comm_output)
      integer, intent(in) :: ndims, rank
      integer, dimension(ndims), intent(in) :: dimensions, processor_dim
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
      integer, dimension(ndims), intent(in) :: processor_dim

      integer, intent(out) :: cart_comm

      call create_cart_communicator_mpi_wrapper(ndims, processor_dim, cart_comm)

   end subroutine create_cart_communicator

   !> Get the cartesian coordinates of a rank in a communicator
   subroutine get_cart_coords(cart_comm, rank, ndims, coords)
      integer, intent(in) :: cart_comm, rank, ndims
      integer, dimension(ndims), intent(inout) :: coords

      call get_cart_coords_mpi_wrapper(cart_comm, rank, ndims, coords)

   end subroutine get_cart_coords

   !> Find the Moore neighborhood of a rank in a cartesian communicator. We should check the order column vs row. Not sure it is correct but it seems to work.
   subroutine get_cart_neighbors(ndims, coords, cart_comm, neighbors, offset_begin, offset_end)
      integer, intent(in) :: ndims, cart_comm
      integer, dimension(ndims), intent(in) :: coords

      integer, dimension(3**ndims), intent(inout) :: neighbors
      integer, dimension(ndims), intent(inout) :: offset_begin, offset_end

      integer :: ii, jj, global_index
      integer, dimension(ndims) :: indices, current_index
      integer :: rank_of_coords, error

      ! Prevent MPI from aborting when we try to find neighbors that do not exist
      call change_MPI_COMM_errhandler_mpi_wrapper(cart_comm)

      ! Determine the ranks of neighboring processes including corners along each dimension.
      offset_begin(:) = 0
      offset_end(:) = 0

      global_index = 1

      if(ndims == 2) then
         do ii = -1,1,1
            do jj = -1,1,1
               current_index = [ii, jj]
               indices = coords + current_index

               call cart_rank_mpi_wrapper(cart_comm, ndims, indices, rank_of_coords, error)

               if(error == 1) then
                  neighbors(global_index) = neighbor_non_existant_rank ! This is a non-existent rank
               else if(ii == 0 .AND. jj == 0) then
                  neighbors(global_index) = neighbor_current_rank ! This is the current rank
               else
                  neighbors(global_index) = rank_of_coords
                  offset_begin(:) = min(offset_begin(:), current_index)
                  offset_end(:) = max(offset_end(:), current_index)
               endif
               global_index = global_index + 1
            end do
         end do
      end if

      ! Restore the original error handler to abort on errors
      call original_MPI_COMM_errhandler_mpi_wrapper(cart_comm)

   end subroutine get_cart_neighbors

   !> Calculate the indices of the elements that need to be sent to the neighbors
   subroutine calculate_neighbor_indices(ndims, dims, send_array, recv_array)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims
      integer, dimension((3**ndims) * 2 * ndims), intent(inout) :: send_array, recv_array

      integer :: ii, jj, global_index
      integer, dimension(2*ndims) :: send_indices, recv_indices

      send_array(:) = -10
      recv_array(:) = -10

      global_index = 1

      if (ndims == 2) then
         do ii = -1,1,1
            call begin_end_send_neighbor_indices(ii, dims(1), send_indices(1), send_indices(ndims + 1))
            call begin_end_recv_neighbor_indices(ii, dims(1), recv_indices(1), recv_indices(ndims + 1))
            do jj = -1,1,1
               call begin_end_send_neighbor_indices(jj, dims(2), send_indices(2), send_indices(ndims + 2))
               call begin_end_recv_neighbor_indices(jj, dims(2), recv_indices(2), recv_indices(ndims + 2))
               send_array(global_index:global_index + 2*ndims-1) = send_indices(:)
               recv_array(global_index:global_index + 2*ndims-1) = recv_indices(:)
               global_index = global_index + 2*ndims
            end do
         end do
      end if

   end subroutine calculate_neighbor_indices

   !> Calculate the indices of the elements that need to be sent to the neighbors from the block
   subroutine begin_end_send_neighbor_indices(index, dim, begin, end)
      integer, intent(in) :: index, dim
      integer, intent(out):: begin, end

      if(index == -1) then
         begin = 2
         end = 2
      else if(index == 0) then
         begin = 2
         end = dim - 1
      else if(index == 1) then
         begin = dim - 1
         end = dim - 1
      endif

   end subroutine begin_end_send_neighbor_indices

   !> Calculate the indices of the elements that need to be written to the block from the neighbors
   subroutine begin_end_recv_neighbor_indices(index, dim, begin, end)
      integer, intent(in) :: index, dim
      integer, intent(out):: begin, end

      if(index == -1) then
         begin = 1
         end = 1
      else if(index == 0) then
         begin = 2
         end = dim - 1
      else if(index == 1) then
         begin = dim
         end = dim
      endif

   end subroutine begin_end_recv_neighbor_indices

   !> Subroutine to get the neighbor indices out
   subroutine get_neighbor_indices(ndims, neighbor_index, indices, begin, end)
      integer, intent(in) :: ndims, neighbor_index
      integer, dimension((3**ndims) * 2 * ndims), intent(in) :: indices
      integer, dimension(ndims), intent(out) :: begin, end

      integer :: start_index

      start_index = (neighbor_index-1) * 2*ndims + 1

      begin(:) = indices(start_index: start_index + ndims - 1)
      end(:) = indices(start_index + ndims: start_index + 2*ndims - 1)

   end subroutine get_neighbor_indices

   !> Subroutine to calculate the number of elements that need to be sent to the neighbors
   subroutine calculate_number_of_elements_to_send(ndims, comm_inout)
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
   subroutine allocate_cart_comm_type(ndims, comm)
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

      integer :: ii, jj, global_index

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "comm_type: "
      write(iounit, *)

      write(iounit, *) "comm: ", comm%comm
      write(iounit, *) "coords: ", comm%coords

      if(ndims == 2) then
         write(iounit, *) "offset_begin: ", comm%offset_begin
         write(iounit, *) "offset_end: ", comm%offset_end
         write(iounit, *) "begin: ", comm%begin
         write(iounit, *) "end: ", comm%end
         write(iounit, *) "size: ", comm%size

         write(iounit, *) "Neighbor element size: "
         do ii = 1, 3
            do jj = 1,3
               global_index = IDX_XD(ndims, [3,3], [ii, jj])
               write(iounit, '(I5)', advance='no') comm%neighbor_elements_size(global_index)
            end do
            write(iounit, *)
         end do

         write(iounit, *) "Neighbor send indices: "
         global_index = 1
         do ii = 1, 3
            do jj = 1, 3
               write(iounit, '(A, I5, I5, I5, I5, A)', advance='no') "[", comm%neighbor_send_indices(global_index), &
                  comm%neighbor_send_indices(global_index + 1), comm%neighbor_send_indices(global_index + 2), &
                  comm%neighbor_send_indices(global_index + 3), "]"
               global_index = global_index + 4
            end do
            write(iounit, *)
         end do

         write(iounit, *) "Neighbor recv indices: "
         global_index = 1
         do ii = 1, 3
            do jj = 1, 3
               write(iounit, '(A, I5, I5, I5, I5, A)', advance='no') "[", comm%neighbor_recv_indices(global_index), &
                  comm%neighbor_recv_indices(global_index + 1), comm%neighbor_recv_indices(global_index + 2), &
                  comm%neighbor_recv_indices(global_index + 3), "]"
               global_index = global_index + 4
            end do
            write(iounit, *)
         end do

         write(iounit, *) "Neighbor ranks: "
         do ii = 1, 3
            do jj = 1, 3
               global_index = IDX_XD(2, [3,3], [ii, jj])
               write(iounit, '(I5)', advance='no') comm%neighbors(global_index)
            end do
            write(iounit, *)
         end do
      endif

      write(iounit, *) "neighbor_send_request: "
      do ii = 1,3
         do jj = 1,3
            global_index = IDX_XD(2, [3,3], [ii, jj])
            write(iounit, '(I5)', advance='no') comm%neighbor_send_request(global_index)
         end do
         write(iounit, *)
      end do

      write(iounit, *) "neighbor_recv_request: "
      do ii = 1,3
         do jj = 1,3
            global_index = IDX_XD(2, [3,3], [ii, jj])
            write(iounit, '(I5)', advance='no') comm%neighbor_recv_request(global_index)
         end do
         write(iounit, *)
      end do

      if(ndims == 3) then
         write(iounit, *) "NOT DONE YET"
      endif

   end subroutine print_cart_comm_type

   !> A routine to print the cartesian grid of the ranks.
   subroutine print_cartesian_grid(cart_comm, world_size, ndim, filename)
      integer, intent(in) :: cart_comm, world_size  ! Cartesian communicator
      integer, intent(in) :: ndim       ! Number of dimensions
      character(len=*), intent(in) :: filename  ! Base filename
      integer :: current_rank, iounit, ios
      integer, dimension(ndim) :: coords  ! Coordinates in the cartesian topology
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
      write(iounit, *) "Rank cartesian processor grid with dimension:", ndim

      ! Iterate over all ranks in the cartesian communicator
      do current_rank = 0, world_size - 1
         ! Get the coordinates for the current rank
         call get_cart_coords_mpi_wrapper(cart_comm, current_rank, ndim, coords)

         ! Write the rank and its coordinates to the file
         write(iounit, '(A, I4, A, *(I4, 1X))') "Rank ", current_rank, " coords: ", coords
      end do

      ! Close the file
      close(iounit)
   end subroutine print_cartesian_grid

end module comm_module
