module print_module
   use rank_parameters_module, only: rank_type
   use utility_functions_module, only: IDX_XD
   use mpi_wrapper_module,
   implicit none

   private
   public :: print_cartesian_grid, print_rank_parameters

contains

   !> A routine to print the cartesian grid of the ranks. For debugging and understanding
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

   !> Print rank parameters
   subroutine print_rank_parameters(parameters, filename)
      type(rank_type), intent(in) :: parameters
      character(255), intent(in) :: filename
      integer :: iounit, ios
      character(255) :: file_with_rank

      ! Create a modified filename by appending the rank to the original filename
      write(file_with_rank, '(A, I0, A)') trim(filename), parameters%rank, ".txt"

      ! Open the file for writing, associate it with a logical unit (iounit)
      open(newunit=iounit, file=file_with_rank, status='replace', action='write', iostat=ios)

      ! Check for errors in opening the file
      if (ios /= 0) then
         print *, "Error opening file: ", file_with_rank
         return
      endif

      ! Write the parameters to the file
      write(iounit, *) "Dimensions:", parameters%ndims
      write(iounit, *) "Rank:", parameters%rank
      write(iounit, *) "World Size:", parameters%world_size
      write(iounit, *) "Grid Size:", parameters%grid_size
      write(iounit, *) "Num Processors per dim:", parameters%processor_dim
      write(iounit, *) "Rank cartesian coordinates:", parameters%coords
      write(iounit, *) "Number of elements to be sent and received:", parameters%num_sendrecv_elements
      write(iounit, *) "Number of elements in the block:", parameters%num_block_elements

      ! Write the block to the file
      call print_block_matrix(parameters%ndims, parameters%block_size, parameters%begin_block, &
         parameters%end_block, parameters%neighbors, parameters%neighbor_sendrecv_start_index, &
         parameters%block_matrix, iounit)

      ! Close the file
      close(iounit)
   end subroutine print_rank_parameters

   !> Print block matrix
   subroutine print_block_matrix(ndims, block_size, begin_block, end_block, block_neighbors, &
      block_neighbor_start_index, block_matrix, iounit)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: block_size, begin_block, end_block
      integer, dimension(3**ndims), intent(in) :: block_neighbors, block_neighbor_start_index
      real, dimension(product(block_size)), intent(in) :: block_matrix
      integer, intent(in) :: iounit

      integer :: ii, jj, kk, global_index

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "Printing block matrix:"
      write(iounit, *)

      write(iounit, *) "Block size:", block_size
      write(iounit, *) "Begin block:", begin_block
      write(iounit, *) "End block:", end_block

      if(ndims == 2) then

         write(iounit, *) "Neighbor ranks:"
         do ii = 1, 3
            do jj = 1, 3
               global_index = IDX_XD(2, [3,3], [jj, ii])
               write(iounit, '(I5)', advance='no') block_neighbors(global_index)
            end do
            write(iounit, *)
         end do

         write(iounit, *) "Neighbor sendrecv start indices:"
         do ii = 1,3
            do jj = 1,3
               global_index = IDX_XD(2, [3,3], [jj, ii])
               write(iounit, '(I5)', advance='no') block_neighbor_start_index(global_index)
            end do
            write(iounit, *)
         end do

         write(iounit, *) "Block matrix:"
         do ii = 1, block_size(1)
            do jj = 1, block_size(2)
               global_index = IDX_XD(2, block_size, [jj, ii])
               write(iounit, '(F10.3)', advance='no') block_matrix(global_index)
            end do
            write(iounit, *)
         end do
      endif

      ! Needs to be updated!
      if(ndims == 3) then
         do ii = 1, block_size(1)
            write(iounit,*) 'Slice (ii = ', ii, '):'
            do jj = 1, block_size(2)
               do kk = 1, block_size(3)
                  global_index = IDX_XD(3, block_size, [kk, jj, ii])
                  write(iounit, '(F10.3)', advance='no') block_matrix(global_index)
               end do
               write(iounit, *)
            end do
            write(iounit, *)
         end do
      endif

   end subroutine print_block_matrix

end module print_module
