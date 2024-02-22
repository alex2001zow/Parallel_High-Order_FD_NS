module block_module
   use constants_module, only: neighbor_non_existant_rank
   use comm_module, only: comm_type
   use neighbor_types_module, only: get_neighbor_range
   use utility_functions_module, only: IDX_XD, sleeper_function
   implicit none

   private

   !> Structure to store the block information.
   type block_type
      integer :: num_elements, num_sendrecv_elements
      integer, allocatable :: size(:), begin(:), end(:)
      real, allocatable :: matrix(:), elements_send(:), elements_recv(:)
      integer, allocatable :: sendrecv_start_index(:)
   end type block_type

   public :: block_type, create_block_type, deallocate_block_type, print_block_type

contains

   !> Subroutine to allocate the block structure.
   subroutine create_block_type(ndims, comm, size, coords, block_output)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: size, coords
      type(comm_type), intent(in) :: comm
      type(block_type), intent(out) :: block_output

      allocate(block_output%size(ndims))
      allocate(block_output%begin(ndims))
      allocate(block_output%end(ndims))

      block_output%size = size

      block_output%begin = coords*size + 1
      block_output%end = block_output%begin + size - 1

      block_output%num_elements = product(size)
      allocate(block_output%matrix(block_output%num_elements))

      block_output%num_sendrecv_elements = 1

      call allocate_neighbor_sendrecv_array(ndims, comm, block_output)

   end subroutine create_block_type

   !> Subroutine to allocate the neighbor sendrecv array.
   subroutine allocate_neighbor_sendrecv_array(ndims, comm, block_inout)
      integer, intent(in) :: ndims
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: block_inout

      integer, dimension(ndims) :: begin, end
      integer :: neighbor_index, neighbor_elements_size

      allocate(block_inout%sendrecv_start_index(comm%num_neighbors))
      block_inout%sendrecv_start_index(:) = -1

      neighbor_elements_size = 1
      do neighbor_index = 1, comm%num_neighbors
         if(comm%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            call get_neighbor_range(ndims, neighbor_index, block_inout%size, begin, end)
            neighbor_elements_size = product(end - begin + 1)
            block_inout%sendrecv_start_index(neighbor_index) = block_inout%num_sendrecv_elements
            block_inout%num_sendrecv_elements = block_inout%num_sendrecv_elements + neighbor_elements_size
         end if
      end do

      block_inout%num_sendrecv_elements = block_inout%num_sendrecv_elements - 1

      allocate(block_inout%elements_send(block_inout%num_sendrecv_elements))
      allocate(block_inout%elements_recv(block_inout%num_sendrecv_elements))

   end subroutine allocate_neighbor_sendrecv_array

   !> Subroutine to deallocate the block structure.
   subroutine deallocate_block_type(block_input)
      type(block_type), intent(inout) :: block_input

      if(allocated(block_input%size)) then
         deallocate(block_input%size)
      endif
      if(allocated(block_input%begin)) then
         deallocate(block_input%begin)
      endif
      if(allocated(block_input%end)) then
         deallocate(block_input%end)
      endif
      if(allocated(block_input%matrix)) then
         deallocate(block_input%matrix)
      endif
      if(allocated(block_input%elements_send)) then
         deallocate(block_input%elements_send)
      endif
      if(allocated(block_input%elements_recv)) then
         deallocate(block_input%elements_recv)
      endif
      if(allocated(block_input%sendrecv_start_index)) then
         deallocate(block_input%sendrecv_start_index)
      endif

   end subroutine deallocate_block_type

   !> Subroutine to print the block structure.
   subroutine print_block_type(ndims, block_input, iounit)
      integer, intent(in) :: ndims, iounit
      type(block_type), intent(in) :: block_input

      integer :: ii, jj, global_index

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "block_type: "
      write(iounit, *)

      write(iounit, *) "size: ", block_input%size
      write(iounit, *) "begin: ", block_input%begin
      write(iounit, *) "end: ", block_input%end
      write(iounit, *) "num_elements: ", block_input%num_elements
      write(iounit, *) "num_sendrecv_elements: ", block_input%num_sendrecv_elements

      if(ndims == 2) then

         write(iounit, *) "sendrecv_start_index: "
         do ii = 1, 3
            do jj = 1, 3
               global_index = IDX_XD(ndims, [3,3], [jj, ii])
               write(iounit, '(I5)', advance='no') block_input%sendrecv_start_index(global_index)
            end do
            write(iounit, *)
         end do

         write(iounit, '(A, *(F10.3))') "elements_send: ", block_input%elements_send(:)
         write(iounit, '(A, *(F10.3))') "elements_recv: ", block_input%elements_recv(:)

         write(iounit, *) "matrix:"
         do ii = 1, block_input%size(1)
            do jj = 1, block_input%size(2)
               global_index = IDX_XD(ndims, block_input%size, [jj, ii])
               write(iounit, '(F10.3)', advance='no') block_input%matrix(global_index)
            end do
            write(iounit, *)
         end do
      endif

      if(ndims == 3) then
         write(iounit, *) "NOT DONE YET"
      endif

   end subroutine print_block_type

end module block_module
