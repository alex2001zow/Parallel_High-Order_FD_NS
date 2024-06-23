module block_module
   use constants_module, only: neighbor_non_existant_rank
   use mpi_wrapper_module, only: block_data_layout_type, define_block_layout, create_send_recv_layout, free_block_layout_type, &
      sendrecv_data_to_neighbors, waitall_mpi_wrapper
   use comm_module, only: comm_type, get_cart_neighbors, begin_end_neighbor_indices
   use utility_functions_module, only: reshape_real_1D_to_2D, reshape_real_1D_to_3D, reshape_real_1D_to_4D, &
      IDX_XD, IDX_XD_INV, calculate_dx, sleeper_function, print_real_1D_array, print_real_2D_array, print_real_3D_array
   implicit none

   private

   !> Structure to store the block information.
   type block_type
      integer :: ndims, elements_per_index, used_elements_per_index
      integer :: num_elements, extended_num_elements
      real, dimension(:), allocatable :: domain_begin, domain_end, grid_dx, extended_grid_dx
      integer, dimension(:), allocatable :: grid_size, extended_grid_size, global_grid_begin, global_grid_end
      integer, dimension(:), allocatable :: ghost_begin, ghost_end
      integer, dimension(:), allocatable :: stencil_begin, stencil_end

      integer, dimension(:), allocatable :: begin_ghost_offset, end_ghost_offset
      integer, dimension(:), allocatable :: begin_stencil_offset, end_stencil_offset

      integer, dimension(:), allocatable :: local_size, extended_local_size

      integer, dimension(:), allocatable :: global_begin_c, global_end_c, global_dims
      integer, dimension(:), allocatable :: extended_global_begin_c, extended_global_end_c, extended_global_dims
      integer, dimension(:), allocatable :: block_begin_c, block_end_c, block_dims
      integer, dimension(:), allocatable :: extended_block_begin_c, extended_block_end_c, extended_block_dims

      ! Arrays for the matrix, f_matrix, and temp_matrix(Jacobi iteration) and the residual_matrix
      real, dimension(:), allocatable :: matrix, f_matrix, temp_matrix, residual_matrix

      ! Direct solver array
      real, dimension(:), allocatable :: direct_solver_matrix
      real, contiguous, dimension(:), pointer :: direct_solver_matrix_ptr
      real, contiguous, dimension(:,:), pointer :: direct_solver_matrix_ptr_2D

      ! The pointers to the matrix, f_matrix, and temp_matrix for different dimensions.
      real, contiguous, dimension(:), pointer :: matrix_ptr, f_matrix_ptr, temp_matrix_ptr, residual_matrix_ptr
      real, contiguous, dimension(:,:), pointer :: matrix_ptr_2D, f_matrix_ptr_2D, temp_matrix_ptr_2D, residual_matrix_ptr_2D
      real, contiguous, dimension(:,:,:), pointer :: matrix_ptr_3D, f_matrix_ptr_3D, temp_matrix_ptr_3D, residual_matrix_ptr_3D
      real, contiguous, dimension(:,:,:,:), pointer :: matrix_ptr_4D, f_matrix_ptr_4D, temp_matrix_ptr_4D, residual_matrix_ptr_4D

      type(block_data_layout_type) :: data_layout
   end type block_type

   public :: block_type, create_block_type, deallocate_block_type, print_block_type, sendrecv_data_neighbors

contains

   !> Subroutine to allocate the block structure. We do not use bc_begin and bc_end yet. Not sure we need them.
   subroutine create_block_type(ndims, elements_per_index, used_elements_per_index, domain_begin, domain_end, grid_size, comm, &
      ghost_begin, ghost_end, stencil_begin, stencil_end, block_output)
      integer, intent(in) :: ndims, elements_per_index, used_elements_per_index
      real, dimension(:), intent(in) :: domain_begin, domain_end
      integer, dimension(:), intent(in) :: grid_size, ghost_begin, ghost_end, stencil_begin, stencil_end
      type(comm_type), intent(in) :: comm
      type(block_type), intent(out) :: block_output

      block_output%ndims = ndims
      block_output%elements_per_index = elements_per_index
      block_output%used_elements_per_index = used_elements_per_index

      call allocate_block_type(block_output)

      block_output%domain_begin = domain_begin
      block_output%domain_end = domain_end
      block_output%grid_size = grid_size
      block_output%ghost_begin = ghost_begin
      block_output%ghost_end = ghost_end
      block_output%stencil_begin = stencil_begin
      block_output%stencil_end = stencil_end
      block_output%extended_grid_size = block_output%grid_size + (block_output%ghost_begin + block_output%ghost_end)

      call calculate_dx(block_output%domain_begin, block_output%domain_end, block_output%grid_size, &
         block_output%grid_dx)

      call calculate_dx(block_output%domain_begin, block_output%domain_end, block_output%extended_grid_size, &
         block_output%extended_grid_dx)

      call setup_block_data_layout_type(block_output%ndims, block_output%elements_per_index, &
         block_output%grid_size, comm%processor_dim, &
         block_output%ghost_begin , block_output%ghost_end, &
         block_output%stencil_begin, block_output%stencil_end, &
         comm%periods, comm%coords, comm%comm, block_output)

      block_output%global_grid_begin = merge(1, 0, comm%periods) * block_output%stencil_begin + block_output%ghost_begin
      block_output%global_grid_end = -1000!block_output%extended_grid_size - &
      !(merge(1, 0, comm%periods) * block_output%stencil_end + block_output%ghost_end)

      allocate(block_output%matrix(block_output%extended_num_elements))
      allocate(block_output%f_matrix(block_output%extended_num_elements))
      allocate(block_output%residual_matrix(block_output%extended_num_elements))

      allocate(block_output%direct_solver_matrix(block_output%extended_num_elements**2))
      block_output%direct_solver_matrix = 0

      block_output%matrix = 0.0
      block_output%f_matrix = 0.0
      block_output%residual_matrix = 0.0


      call setup_matrix_pointers(block_output)

   end subroutine create_block_type

   ! Setup the block data layout type using 0-based indexing due to OpenMPI.
   subroutine setup_block_data_layout_type(ndims, elements_per_index, grid_size, processor_dims, &
      ghost_begin, ghost_end, stencil_begin, stencil_end, &
      periods, coords, comm, block_inout)
      integer, intent(in) :: ndims, elements_per_index, comm
      integer, dimension(:), intent(in) :: grid_size, processor_dims
      integer, dimension(:), intent(in) :: ghost_begin, ghost_end, stencil_begin, stencil_end, coords
      logical, dimension(:), intent(in) :: periods
      type(block_type), target, intent(inout) :: block_inout

      type(block_data_layout_type), pointer:: block_data_layout

      integer, dimension(ndims) :: begin_ghost_offset, end_ghost_offset
      integer, dimension(ndims) :: begin_stencil_offset, end_stencil_offset

      integer, dimension(ndims) :: local_size, extended_local_size

      integer, dimension(ndims) :: global_begin_c, global_end_c, global_dims
      integer, dimension(ndims) :: extended_global_begin_c, extended_global_end_c, extended_global_dims
      integer, dimension(ndims) :: block_begin_c, block_end_c, block_dims
      integer, dimension(ndims) :: extended_block_begin_c, extended_block_end_c, extended_block_dims

      integer, dimension(3**ndims) :: neighbors
      integer, dimension(ndims) :: begin_send, end_send, send_dims
      integer, dimension(ndims) :: begin_recv, end_recv, recv_dims
      integer :: true_neighbor_index, global_index
      integer, dimension(ndims) :: stencil_size, current_index, neighbor_array
      integer :: recv_tag

      block_data_layout => block_inout%data_layout

      !> Find the neighbors of the current block
      call get_cart_neighbors(ndims, coords, comm, neighbors, &
         begin_ghost_offset, end_ghost_offset, begin_stencil_offset, end_stencil_offset)

      ! Allocate the neighbor ranks and send/recv types
      block_data_layout%number_of_existing_neighbors = count(neighbors > neighbor_non_existant_rank)
      allocate(block_data_layout%neighbor_ranks(block_data_layout%number_of_existing_neighbors))
      allocate(block_data_layout%send_type_4(block_data_layout%number_of_existing_neighbors))
      allocate(block_data_layout%recv_type_4(block_data_layout%number_of_existing_neighbors))
      allocate(block_data_layout%send_tag_4(block_data_layout%number_of_existing_neighbors))
      allocate(block_data_layout%recv_tag_4(block_data_layout%number_of_existing_neighbors))

      ! Scale the ghost by the ghost sizes
      begin_ghost_offset = begin_ghost_offset * ghost_begin
      end_ghost_offset = end_ghost_offset * ghost_end

      ! Scale the stencil by the stencil sizes
      begin_stencil_offset = begin_stencil_offset * stencil_begin
      end_stencil_offset = end_stencil_offset * stencil_end

      ! Calculate the local size and extended local size
      local_size = grid_size / processor_dims
      extended_local_size = local_size + (begin_ghost_offset + end_ghost_offset) + (begin_stencil_offset + end_stencil_offset)

      ! Calculate the global begin and end coordinates using 0-based indexing due to OpenMPI
      global_begin_c = coords * local_size + (begin_ghost_offset + begin_stencil_offset)
      global_end_c = global_begin_c + local_size
      global_dims = global_end_c - global_begin_c

      ! Calculate the extended global begin and end coordinates. Not sure if correct
      extended_global_begin_c = global_begin_c - (begin_ghost_offset + begin_stencil_offset)
      extended_global_end_c = global_end_c + (end_ghost_offset + end_stencil_offset)
      extended_global_dims = extended_global_end_c - extended_global_begin_c

      ! Calculate the block begin and end coordinates
      block_begin_c = begin_ghost_offset + begin_stencil_offset
      block_end_c = block_begin_c + local_size
      block_dims = block_end_c - block_begin_c

      ! Calculate the extended block begin and end coordinates
      extended_block_begin_c = 0
      extended_block_end_c = extended_block_begin_c + extended_local_size
      extended_block_dims = extended_block_end_c - extended_block_begin_c

      ! Calculate the number of elements and extended number of elements
      block_inout%num_elements = product(block_end_c - block_begin_c) * elements_per_index
      block_inout%extended_num_elements = product(extended_block_end_c - extended_block_begin_c) * elements_per_index

      ! Set the values in the block
      block_inout%begin_ghost_offset = begin_ghost_offset
      block_inout%end_ghost_offset = end_ghost_offset
      block_inout%begin_stencil_offset = begin_stencil_offset
      block_inout%end_stencil_offset = end_stencil_offset

      block_inout%local_size = local_size
      block_inout%extended_local_size = extended_local_size

      block_inout%global_begin_c = global_begin_c
      block_inout%global_end_c = global_end_c
      block_inout%global_dims = global_dims

      block_inout%extended_global_begin_c = extended_global_begin_c
      block_inout%extended_global_end_c = extended_global_end_c
      block_inout%extended_global_dims = extended_global_dims

      block_inout%block_begin_c = block_begin_c
      block_inout%block_end_c = block_end_c
      block_inout%block_dims = block_dims

      block_inout%extended_block_begin_c = extended_block_begin_c
      block_inout%extended_block_end_c = extended_block_end_c
      block_inout%extended_block_dims = extended_block_dims

      ! Define the block data layout
      call define_block_layout(ndims, elements_per_index, block_inout%grid_size, &
         ghost_begin + merge(1, 0, periods) * stencil_begin, &
         ghost_end + merge(1, 0, periods) * stencil_end, &
         global_begin_c, global_dims, &
         extended_global_begin_c, extended_global_dims, &
         block_begin_c, block_dims, &
         extended_block_begin_c, extended_block_dims, &
         block_data_layout)

      stencil_size = (stencil_begin + stencil_end)/2

      ! Define the send and recv types for the neighbors
      neighbor_array = 3 ! 3 is the number of neighbors per dimension

      true_neighbor_index = 1
      do global_index = 1, 3**ndims
         call IDX_XD_INV(ndims, neighbor_array, global_index, current_index)
         current_index = current_index - 2 ! To go -1, 0, 1

         if(neighbors(global_index) > neighbor_non_existant_rank) then
            block_data_layout%neighbor_ranks(true_neighbor_index) = neighbors(global_index)

            ! Send part
            call begin_end_neighbor_indices(ndims, extended_block_begin_c, extended_block_end_c, &
               block_begin_c, block_end_c, stencil_size, current_index, 0, begin_send, end_send, send_dims)
            block_data_layout%send_tag_4(true_neighbor_index) = global_index

            ! Recv part
            call begin_end_neighbor_indices(ndims, extended_block_begin_c, extended_block_end_c, &
               block_begin_c, block_end_c, stencil_size, current_index, 1, begin_recv, end_recv, recv_dims)
            call IDX_XD(ndims,neighbor_array, (-current_index+2), recv_tag)
            block_data_layout%recv_tag_4(true_neighbor_index) = recv_tag

            !print *, (current_index+2), (-current_index+2), global_index, recv_tag

            call create_send_recv_layout(ndims, true_neighbor_index, extended_local_size, &
               begin_send, end_send, begin_recv, end_recv, block_data_layout)

            true_neighbor_index = true_neighbor_index + 1
         end if
      end do

   end subroutine setup_block_data_layout_type

   !> Subroutine to send and recieve between neighbors.
   subroutine sendrecv_data_neighbors(comm, block_input, matrix)
      integer, intent(in) :: comm
      type(block_type), intent(inout) :: block_input
      real, dimension(:), intent(inout) :: matrix

      integer :: neighbor_index, neighbor_rank
      integer, dimension(block_input%data_layout%number_of_existing_neighbors*2) :: sendrecv_request

      do neighbor_index = 1, block_input%data_layout%number_of_existing_neighbors
         neighbor_rank = block_input%data_layout%neighbor_ranks(neighbor_index)
         call sendrecv_data_to_neighbors(block_input%data_layout, matrix, comm, neighbor_rank, neighbor_index, &
            sendrecv_request(neighbor_index), &
            sendrecv_request(neighbor_index + block_input%data_layout%number_of_existing_neighbors))
      end do

      call waitall_mpi_wrapper(block_input%data_layout%number_of_existing_neighbors*2, sendrecv_request)

   end subroutine sendrecv_data_neighbors

   !> Setup pointers to the matrix, f_matrix, and temp_matrix and create reshaped pointers for 2D, 3D, and 4D matrices. Can be extended to 5D, 6D, etc.
   subroutine setup_matrix_pointers(block_input)
      type(block_type), intent(inout), target :: block_input

      ! Set the pointers to the matrix, f_matrix, and temp_matrix
      block_input%matrix_ptr => block_input%matrix
      block_input%f_matrix_ptr => block_input%f_matrix
      block_input%temp_matrix_ptr => block_input%temp_matrix
      block_input%residual_matrix_ptr => block_input%residual_matrix

      ! Set the pointers to the direct solver matrix
      block_input%direct_solver_matrix_ptr => block_input%direct_solver_matrix

      ! Reshape the pointers to the matrix, f_matrix, and temp_matrix for 2D, 3D, and 4D matrices
      if(block_input%ndims == 2) then
         call reshape_real_1D_to_2D(block_input%extended_block_dims, block_input%matrix_ptr, &
            block_input%matrix_ptr_2D)
         call reshape_real_1D_to_2D(block_input%extended_block_dims, block_input%f_matrix_ptr, &
            block_input%f_matrix_ptr_2D)
         call reshape_real_1D_to_2D(block_input%extended_block_dims, block_input%temp_matrix_ptr, &
            block_input%temp_matrix_ptr_2D)
         call reshape_real_1D_to_2D(block_input%extended_block_dims, block_input%residual_matrix_ptr, &
            block_input%residual_matrix_ptr_2D)

         !!! TEST
         call reshape_real_1D_to_2D([block_input%extended_block_dims(1)*block_input%extended_block_dims(2), &
            block_input%extended_block_dims(2)*block_input%extended_block_dims(1)], block_input%direct_solver_matrix_ptr, &
            block_input%direct_solver_matrix_ptr_2D)

      else if(block_input%ndims == 3) then
         call reshape_real_1D_to_3D(block_input%extended_block_dims, block_input%matrix_ptr, &
            block_input%matrix_ptr_3D)
         call reshape_real_1D_to_3D(block_input%extended_block_dims, block_input%f_matrix_ptr, &
            block_input%f_matrix_ptr_3D)
         call reshape_real_1D_to_3D(block_input%extended_block_dims, block_input%temp_matrix_ptr, &
            block_input%temp_matrix_ptr_3D)
         call reshape_real_1D_to_3D(block_input%extended_block_dims, block_input%residual_matrix_ptr, &
            block_input%residual_matrix_ptr_3D)
      else if(block_input%ndims == 4) then
         call reshape_real_1D_to_4D(block_input%extended_block_dims, block_input%matrix_ptr, &
            block_input%matrix_ptr_4D)
         call reshape_real_1D_to_4D(block_input%extended_block_dims, block_input%f_matrix_ptr, &
            block_input%f_matrix_ptr_4D)
         call reshape_real_1D_to_4D(block_input%extended_block_dims, block_input%temp_matrix_ptr, &
            block_input%temp_matrix_ptr_4D)
         call reshape_real_1D_to_4D(block_input%extended_block_dims, block_input%residual_matrix_ptr, &
            block_input%residual_matrix_ptr_4D)
      end if

   end subroutine setup_matrix_pointers

   !> Subroutine to allocate the block structure.
   pure subroutine allocate_block_type(block_input)
      type(block_type), intent(inout) :: block_input

      allocate(block_input%domain_begin(block_input%ndims))
      allocate(block_input%domain_end(block_input%ndims))
      allocate(block_input%grid_dx(block_input%ndims))
      allocate(block_input%extended_grid_dx(block_input%ndims))
      allocate(block_input%grid_size(block_input%ndims))
      allocate(block_input%extended_grid_size(block_input%ndims))
      allocate(block_input%global_grid_begin(block_input%ndims))
      allocate(block_input%global_grid_end(block_input%ndims))
      allocate(block_input%ghost_begin(block_input%ndims))
      allocate(block_input%ghost_end(block_input%ndims))
      allocate(block_input%stencil_begin(block_input%ndims))
      allocate(block_input%stencil_end(block_input%ndims))

      allocate(block_input%begin_ghost_offset(block_input%ndims))
      allocate(block_input%end_ghost_offset(block_input%ndims))
      allocate(block_input%begin_stencil_offset(block_input%ndims))
      allocate(block_input%end_stencil_offset(block_input%ndims))

      allocate(block_input%local_size(block_input%ndims))
      allocate(block_input%extended_local_size(block_input%ndims))

      allocate(block_input%global_begin_c(block_input%ndims))
      allocate(block_input%global_end_c(block_input%ndims))
      allocate(block_input%global_dims(block_input%ndims))

      allocate(block_input%extended_global_begin_c(block_input%ndims))
      allocate(block_input%extended_global_end_c(block_input%ndims))
      allocate(block_input%extended_global_dims(block_input%ndims))

      allocate(block_input%block_begin_c(block_input%ndims))
      allocate(block_input%block_end_c(block_input%ndims))
      allocate(block_input%block_dims(block_input%ndims))

      allocate(block_input%extended_block_begin_c(block_input%ndims))
      allocate(block_input%extended_block_end_c(block_input%ndims))
      allocate(block_input%extended_block_dims(block_input%ndims))

   end subroutine allocate_block_type

   !> Subroutine to deallocate the block structure.
   subroutine deallocate_block_type(block_input)
      type(block_type), intent(inout) :: block_input

      if(allocated(block_input%domain_begin)) deallocate(block_input%domain_begin)
      if(allocated(block_input%domain_end)) deallocate(block_input%domain_end)
      if(allocated(block_input%grid_dx)) deallocate(block_input%grid_dx)
      if(allocated(block_input%extended_grid_dx)) deallocate(block_input%extended_grid_dx)
      if(allocated(block_input%grid_size)) deallocate(block_input%grid_size)
      if(allocated(block_input%extended_grid_size)) deallocate(block_input%extended_grid_size)
      if(allocated(block_input%global_grid_begin)) deallocate(block_input%global_grid_begin)
      if(allocated(block_input%global_grid_end)) deallocate(block_input%global_grid_end)
      if(allocated(block_input%ghost_begin)) deallocate(block_input%ghost_begin)
      if(allocated(block_input%ghost_end)) deallocate(block_input%ghost_end)
      if(allocated(block_input%stencil_begin)) deallocate(block_input%stencil_begin)
      if(allocated(block_input%stencil_end)) deallocate(block_input%stencil_end)

      if(allocated(block_input%begin_ghost_offset)) deallocate(block_input%begin_ghost_offset)
      if(allocated(block_input%end_ghost_offset)) deallocate(block_input%end_ghost_offset)
      if(allocated(block_input%begin_stencil_offset)) deallocate(block_input%begin_stencil_offset)
      if(allocated(block_input%end_stencil_offset)) deallocate(block_input%end_stencil_offset)

      if(allocated(block_input%local_size)) deallocate(block_input%local_size)
      if(allocated(block_input%extended_local_size)) deallocate(block_input%extended_local_size)

      if(allocated(block_input%global_begin_c)) deallocate(block_input%global_begin_c)
      if(allocated(block_input%global_end_c)) deallocate(block_input%global_end_c)
      if(allocated(block_input%global_dims)) deallocate(block_input%global_dims)

      if(allocated(block_input%extended_global_begin_c)) deallocate(block_input%extended_global_begin_c)
      if(allocated(block_input%extended_global_end_c)) deallocate(block_input%extended_global_end_c)
      if(allocated(block_input%extended_global_dims)) deallocate(block_input%extended_global_dims)

      if(allocated(block_input%block_begin_c)) deallocate(block_input%block_begin_c)
      if(allocated(block_input%block_end_c)) deallocate(block_input%block_end_c)
      if(allocated(block_input%block_dims)) deallocate(block_input%block_dims)

      if(allocated(block_input%extended_block_begin_c)) deallocate(block_input%extended_block_begin_c)
      if(allocated(block_input%extended_block_end_c)) deallocate(block_input%extended_block_end_c)
      if(allocated(block_input%extended_block_dims)) deallocate(block_input%extended_block_dims)

      if(allocated(block_input%matrix)) deallocate(block_input%matrix)
      if(allocated(block_input%f_matrix)) deallocate(block_input%f_matrix)
      if(allocated(block_input%temp_matrix)) deallocate(block_input%temp_matrix)
      if(allocated(block_input%residual_matrix)) deallocate(block_input%residual_matrix)

      !!! TEST
      if(allocated(block_input%direct_solver_matrix)) deallocate(block_input%direct_solver_matrix)

      call free_block_layout_type(block_input%data_layout)

   end subroutine deallocate_block_type

   !> Subroutine to print the block structure.
   subroutine print_block_type(block_input, iounit)
      type(block_type), intent(in) :: block_input
      integer, intent(in) :: iounit

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "block_type: "
      write(iounit, *)

      write(iounit, *) "ndims: ", block_input%ndims
      write(iounit, *) "elements_per_index: ", block_input%elements_per_index
      write(iounit, *) "used_elements_per_index: ", block_input%used_elements_per_index

      write(iounit, *)

      write(iounit, *) "domain_begin: ", block_input%domain_begin
      write(iounit, *) "domain_end: ", block_input%domain_end
      write(iounit, *) "grid_dx: ", block_input%grid_dx
      write(iounit, *) "extended_grid_dx: ", block_input%extended_grid_dx
      write(iounit, *) "grid_size: ", block_input%grid_size
      write(iounit, *) "extended_grid_size: ", block_input%extended_grid_size
      write(iounit, *) "global_grid_begin: ", block_input%global_grid_begin
      write(iounit, *) "global_grid_end: ", block_input%global_grid_end
      write(iounit, *) "ghost_begin: ", block_input%ghost_begin
      write(iounit, *) "ghost_end: ", block_input%ghost_end
      write(iounit, *) "stencil_begin: ", block_input%stencil_begin
      write(iounit, *) "stencil_end: ", block_input%stencil_end

      write(iounit, *)

      write(iounit, *) "begin_ghost_offset: ", block_input%begin_ghost_offset
      write(iounit, *) "end_ghost_offset: ", block_input%end_ghost_offset
      write(iounit, *) "begin_stencil_offset: ", block_input%begin_stencil_offset
      write(iounit, *) "end_stencil_offset: ", block_input%end_stencil_offset

      write(iounit, *)

      write(iounit, *) "local_size: ", block_input%local_size
      write(iounit, *) "extended_local_size: ", block_input%extended_local_size

      write(iounit, *)

      write(iounit, *) "global_begin_c: ", block_input%global_begin_c
      write(iounit, *) "global_end_c: ", block_input%global_end_c
      write(iounit, *) "global_dims: ", block_input%global_dims

      write(iounit, *)

      write(iounit, *) "extended_global_begin_c: ", block_input%extended_global_begin_c
      write(iounit, *) "extended_global_end_c: ", block_input%extended_global_end_c
      write(iounit, *) "extended_global_dims: ", block_input%extended_global_dims

      write(iounit, *)

      write(iounit, *) "block_begin_c: ", block_input%block_begin_c
      write(iounit, *) "block_end_c: ", block_input%block_end_c
      write(iounit, *) "block_dims: ", block_input%block_dims

      write(iounit, *)

      write(iounit, *) "extended_block_begin_c: ", block_input%extended_block_begin_c
      write(iounit, *) "extended_block_end_c: ", block_input%extended_block_end_c
      write(iounit, *) "extended_block_dims: ", block_input%extended_block_dims

      write(iounit, *)

      write(iounit, *) "num_elements: ", block_input%num_elements
      write(iounit, *) "extended_num_elements: ", block_input%extended_num_elements

      write(iounit, *)
      write(iounit, *) "Matrix"
      select case(block_input%ndims)
       case (1)
         call print_real_1D_array(block_input%matrix_ptr, iounit)
       case (2)
         call print_real_2D_array(block_input%matrix_ptr_2D, iounit)
       case (3)
         call print_real_3D_array(block_input%matrix_ptr_3D, iounit)
      end select

      write(iounit, *)
      write(iounit, *) "F Matrix"
      select case(block_input%ndims)
       case (1)
         call print_real_1D_array(block_input%f_matrix_ptr, iounit)
       case (2)
         call print_real_2D_array(block_input%f_matrix_ptr_2D, iounit)
       case (3)
         call print_real_3D_array(block_input%f_matrix_ptr_3D, iounit)
      end select

      ! write(iounit, *)
      ! write(iounit, *) "Residual Matrix"
      ! select case(block_input%ndims)
      !  case (1)
      !    call print_real_1D_array(block_input%residual_matrix_ptr, iounit)
      !  case (2)
      !    call print_real_2D_array(block_input%residual_matrix_ptr_2D, iounit)
      !  case (3)
      !    call print_real_3D_array(block_input%residual_matrix_ptr_3D, iounit)
      ! end select
      !call print_real_2D_array(block_input%direct_solver_matrix_ptr_2D, iounit)

   end subroutine print_block_type

end module block_module
