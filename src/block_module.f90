module block_module
   use constants_module, only: neighbor_non_existant_rank
   use comm_module, only: comm_type
   use utility_functions_module, only: IDX_XD, IDX_XD_INV, print_real_array, print_integer_array, sleeper_function
   implicit none

   private

   !> Structure to store the block information.
   type block_type
      integer :: num_elements, num_sendrecv_elements
      integer, dimension(:), allocatable :: size, begin, end
      real, dimension(:), allocatable :: matrix, f_array, temp_array
      real, dimension(:), allocatable :: elements_send, elements_recv
      integer, dimension(:), allocatable :: sendrecv_start_index
   end type block_type

   public :: block_type, create_block_type, deallocate_block_type, print_block_type

contains

   !> Subroutine to allocate the block structure.
   pure subroutine create_block_type(ndims, comm, block_output)
      integer, intent(in) :: ndims
      type(comm_type), intent(in) :: comm
      type(block_type), intent(out) :: block_output

      allocate(block_output%size(ndims))
      allocate(block_output%begin(ndims))
      allocate(block_output%end(ndims))

      block_output%begin = comm%begin
      block_output%end = comm%end

      block_output%size = comm%size

      block_output%num_elements = product(block_output%size)
      allocate(block_output%matrix(block_output%num_elements))

      call allocate_neighbor_sendrecv_array(comm, block_output)

   end subroutine create_block_type

   !> Subroutine to allocate the neighbor sendrecv array.
   pure subroutine allocate_neighbor_sendrecv_array(comm, block_inout)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: block_inout

      integer :: neighbor_index, neighbor_elements_size

      block_inout%num_sendrecv_elements = 1

      allocate(block_inout%sendrecv_start_index(comm%num_neighbors))
      block_inout%sendrecv_start_index(:) = -1

      neighbor_elements_size = 1
      do neighbor_index = 1, comm%num_neighbors
         if(comm%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            neighbor_elements_size = comm%neighbor_elements_size(neighbor_index)
            block_inout%sendrecv_start_index(neighbor_index) = block_inout%num_sendrecv_elements
            block_inout%num_sendrecv_elements = block_inout%num_sendrecv_elements + neighbor_elements_size
         end if
      end do

      block_inout%num_sendrecv_elements = block_inout%num_sendrecv_elements - 1

      allocate(block_inout%elements_send(block_inout%num_sendrecv_elements))
      allocate(block_inout%elements_recv(block_inout%num_sendrecv_elements))

      block_inout%elements_send(:) = 10.0
      block_inout%elements_recv(:) = 20.0

   end subroutine allocate_neighbor_sendrecv_array

   !> Subroutine to deallocate the block structure.
   pure subroutine deallocate_block_type(block_input)
      type(block_type), intent(inout) :: block_input

      if(allocated(block_input%size)) deallocate(block_input%size)
      if(allocated(block_input%begin)) deallocate(block_input%begin)
      if(allocated(block_input%end)) deallocate(block_input%end)

      if(allocated(block_input%matrix)) deallocate(block_input%matrix)
      if(allocated(block_input%elements_send)) deallocate(block_input%elements_send)
      if(allocated(block_input%elements_recv)) deallocate(block_input%elements_recv)
      if(allocated(block_input%f_array)) deallocate(block_input%f_array)
      if(allocated(block_input%temp_array)) deallocate(block_input%temp_array)

      if(allocated(block_input%sendrecv_start_index)) deallocate(block_input%sendrecv_start_index)

   end subroutine deallocate_block_type

   !> Subroutine to print the block structure.
   subroutine print_block_type(ndims, block_input, iounit)
      integer, intent(in) :: ndims, iounit
      type(block_type), intent(in) :: block_input

      integer, dimension(ndims) :: neighbor_array

      neighbor_array(:) = 3 ! 3 is the number of neighbors per dimension

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "block_type: "
      write(iounit, *)

      write(iounit, *) "size: ", block_input%size
      write(iounit, *) "begin: ", block_input%begin
      write(iounit, *) "end: ", block_input%end
      write(iounit, *) "num_elements: ", block_input%num_elements
      write(iounit, *) "num_sendrecv_elements: ", block_input%num_sendrecv_elements

      call print_integer_array(ndims, neighbor_array, block_input%sendrecv_start_index, 1, "sendrecv_start_index: ", iounit)

      write(iounit, '(A, *(F10.3))') "elements_send: ", block_input%elements_send(:)
      write(iounit, '(A, *(F10.3))') "elements_recv: ", block_input%elements_recv(:)

      call print_real_array(ndims, block_input%size, block_input%matrix, 1, "Matrix", iounit)

   end subroutine print_block_type

end module block_module
