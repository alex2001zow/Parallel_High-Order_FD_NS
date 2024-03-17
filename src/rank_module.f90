module rank_module
   use constants_module, only: neighbor_current_rank, neighbor_non_existant_rank
   use functions_module, only: FunctionPair, set_function_pointers
   use mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_PRECISION, MPI_INT
   use comm_module, only : comm_type, create_cart_comm_type, deallocate_cart_comm_type, print_cart_comm_type, &
      get_neighbor_indices
   use block_module, only : block_type, create_block_type, deallocate_block_type, print_block_type
   use finite_difference_module, only : FDstencil_type, create_finite_difference_stencil, &
      deallocate_finite_difference_stencil, print_finite_difference_stencil
   use mpi_wrapper_module, only: create_cart_communicator_mpi_wrapper, get_cart_coords_mpi_wrapper, &
      cart_rank_mpi_wrapper, change_MPI_COMM_errhandler_mpi_wrapper, &
      original_MPI_COMM_errhandler_mpi_wrapper, isendrecv_mpi_wrapper, waitall_mpi_wrapper, free_communicator_mpi_wrapper, &
      write_to_file_mpi_wrapper
   use utility_functions_module, only : IDX_XD, sleeper_function
   implicit none

   private

   !> Structure to hold the parameters for each rank.
   type rank_type
      integer :: ndims, rank, world_size
      integer, dimension(:), allocatable :: grid_size, processor_dim
      real, dimension(:), allocatable :: domain_begin, domain_end

      type(FunctionPair) :: funcs
      type(comm_type) :: comm
      type(FDstencil_type) :: FDstencil
      type(block_type) :: block

   end type rank_type

   public :: rank_type
   public :: create_rank_type, deallocate_rank_type
   public :: communicate_step
   public :: print_rank_type, write_rank_type_blocks_to_file

contains

   !> Routine to create the rank_type
   subroutine create_rank_type(ndims, grid_size, processor_dim, func_enum, rank, world_size, parameters)
      integer, intent(in) :: ndims, rank, world_size, func_enum
      integer, dimension(:), intent(in) :: grid_size, processor_dim
      type(rank_type), intent(inout) :: parameters

      integer, parameter :: num_derivatives = 2
      integer, dimension(2*num_derivatives), parameter :: derivatives = [1,3,3,1]
      integer, dimension(2*num_derivatives), parameter :: derivatives_order = [0,2,2,0]
      real, dimension(num_derivatives), parameter :: derivatives_sign = [1,1]

      integer, dimension(2), parameter :: alphas = [1,1], betas = [1,1]
      real, dimension(2) :: dx

      parameters%ndims = ndims
      parameters%rank = rank
      parameters%world_size = world_size

      call allocate_rank_type(parameters)
      parameters%grid_size(:) = grid_size(:)
      parameters%processor_dim(:) = processor_dim(:)

      parameters%domain_begin(:) = 0.0
      parameters%domain_end(:) = 1.0

      dx = abs((parameters%domain_end(:) - parameters%domain_begin(:))) / (parameters%grid_size(:)-1)

      call set_function_pointers(func_enum, parameters%funcs)

      call create_cart_comm_type(parameters%ndims, parameters%grid_size, parameters%processor_dim, parameters%rank, parameters%comm)

      call create_finite_difference_stencil(parameters%ndims, num_derivatives, derivatives, derivatives_order,&
         derivatives_sign, dx, alphas, betas, parameters%FDstencil)

      call create_block_type(parameters%ndims, parameters%comm, parameters%block)

   end subroutine create_rank_type

   !> Routine to allocate the rank_type
   subroutine allocate_rank_type(parameters)
      type(rank_type), intent(inout) :: parameters

      allocate(parameters%grid_size(parameters%ndims))
      allocate(parameters%processor_dim(parameters%ndims))
      allocate(parameters%domain_begin(parameters%ndims))
      allocate(parameters%domain_end(parameters%ndims))

   end subroutine allocate_rank_type

   !> Routine to deallocate the rank_type
   subroutine deallocate_rank_type(parameters)
      type(rank_type), intent(inout) :: parameters

      if (allocated(parameters%grid_size)) deallocate(parameters%grid_size)
      if (allocated(parameters%processor_dim)) deallocate(parameters%processor_dim)
      if (allocated(parameters%domain_begin)) deallocate(parameters%domain_begin)
      if (allocated(parameters%domain_end)) deallocate(parameters%domain_end)

      call deallocate_cart_comm_type(parameters%comm)
      call deallocate_finite_difference_stencil(parameters%FDstencil)
      call deallocate_block_type(parameters%block)

   end subroutine deallocate_rank_type

   !> Routine to print the rank_type
   subroutine print_rank_type(parameters, filename)
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
      end if

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "rank_type: "
      write(iounit, *)

      write(iounit, *) "ndims:", parameters%ndims
      write(iounit, *) "rank: ", parameters%rank
      write(iounit, *) "world_size: ", parameters%world_size

      write(iounit, *) "grid_size: ", parameters%grid_size
      write(iounit, *) "processor_dim: ", parameters%processor_dim

      call print_cart_comm_type(parameters%ndims, parameters%comm, iounit)

      call print_finite_difference_stencil(parameters%ndims, parameters%FDstencil, iounit)

      call print_block_type(parameters%ndims, parameters%block, iounit)

      ! Close the file
      close(iounit)

   end subroutine print_rank_type

   !> Routine to write the block data to a file
   subroutine write_rank_type_blocks_to_file(parameters, filename_system_solution)
      character(255), intent(in) :: filename_system_solution
      type(rank_type), intent(in) :: parameters

      integer, dimension(parameters%ndims) :: block_dims, block_begin
      integer :: offset_solution, elements_to_write, starting_offset, extent
      integer :: num_blocks_type_vector, whole_block_count

      integer :: size_real

      size_real = sizeof(0.0d0) ! Size of a double precision number in bytes

      block_dims = parameters%block%size ! Size of the block (including ghost points)
      block_begin = parameters%block%begin ! Start index of the block (including ghost points)

      num_blocks_type_vector = 1
      whole_block_count = product(block_dims)

      elements_to_write = whole_block_count ! Number of elements to write to the file
      offset_solution = (parameters%rank) * elements_to_write * size_real ! Where in the file to start writing the true block. Times 8 for double precision
      starting_offset = product(block_begin) * size_real ! The starting offset for the true block in the full block. This is for MPI_TYPE_VECTOR.
      extent = elements_to_write * size_real ! Where the true block ends in the full block. Times 8 for double precision

      call write_to_file_mpi_wrapper(parameters%ndims, & ! number of dimensions
         filename_system_solution, & ! name of the file
         parameters%comm%comm, & ! MPI communicator
         num_blocks_type_vector, & ! num blocks for MPI_TYPE_VECTOR. This is the total number of blocks (including ghost points)
         whole_block_count, & ! block_length (no ghost points)
         whole_block_count, & ! full block length (stride) (including ghost points)
         0, & ! starting offset for MPI_TYPE_RESIZED. Currently found using MPI_GET_EXTENT for just the whole block.
         0, & ! extent for the MPI_TYPE_RESIZED. Currently found using MPI_GET_EXTENT for just the whole block.
         parameters%grid_size, & ! global dims
         block_dims, & ! local dims of the true block (no ghost points)
         block_begin - 1, & ! start index of the true block (no ghost points). Minus 1 to convert to 0 based indexing for MPI
         parameters%block%matrix, & ! the block matrix
         0, & ! offset in the file to start writing the true block
         elements_to_write) ! number of elements to write to the file
   end subroutine write_rank_type_blocks_to_file

   !> Routine to communicate with the neighbors
   !! We could edit this routine to initate the recieve first then write to the send buffer, initate send and then check if recv then write to buffer then wait for send and done.
   subroutine communicate_step(ndims, comm, block, matrix)
      integer, intent(in) :: ndims
      type(comm_type), intent(inout) :: comm
      type(block_type), intent(inout) :: block
      real, dimension(product(block%size)), intent(inout) :: matrix

      integer, dimension(ndims) :: begin, end
      integer :: neighbor_index

      ! Write the data to the send buffer
      do neighbor_index = 1, comm%num_neighbors
         if(comm%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            call get_neighbor_indices(ndims, neighbor_index, comm%neighbor_send_indices, begin, end)
            call array2buffer(ndims, block%size, begin, end, matrix, block%num_elements, block%elements_send, &
               block%num_sendrecv_elements, block%sendrecv_start_index(neighbor_index))
         end if
      end do

      ! Send and recieve the data from the neighbors
      call sendrecv_elements_from_neighbors(comm, block)

      ! Wait for all the sends and recieves to complete
      call waitall_mpi_wrapper(comm%num_neighbors, comm%neighbor_send_request)
      call waitall_mpi_wrapper(comm%num_neighbors, comm%neighbor_recv_request)

      ! Write the data from the recv buffer to the array
      do neighbor_index = 1, comm%num_neighbors
         if(comm%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            call get_neighbor_indices(ndims,neighbor_index, comm%neighbor_recv_indices, begin, end)
            call buffer2array(ndims, block%size, begin, end, matrix, block%num_elements, block%elements_recv, &
               block%num_sendrecv_elements, block%sendrecv_start_index(neighbor_index))
         end if
      end do

   end subroutine communicate_step

   !> Routine to send and recieve elements from each neighbor
   subroutine sendrecv_elements_from_neighbors(comm, block)
      type(comm_type), intent(inout) :: comm
      type(block_type), intent(inout) :: block

      integer :: neighbor_index, neighbor_elements_size

      do neighbor_index = 1, comm%num_neighbors
         if(comm%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            neighbor_elements_size = comm%neighbor_elements_size(neighbor_index)
            call isendrecv_mpi_wrapper(block%num_sendrecv_elements, block%elements_send, &
               block%elements_recv, block%sendrecv_start_index(neighbor_index), &
               neighbor_elements_size, comm%neighbors(neighbor_index),&
               comm%comm, comm%neighbor_send_request(neighbor_index), &
               comm%neighbor_recv_request(neighbor_index))
         else
            comm%neighbor_send_request(neighbor_index) = MPI_REQUEST_NULL
            comm%neighbor_recv_request(neighbor_index) = MPI_REQUEST_NULL
         end if
      end do

   end subroutine sendrecv_elements_from_neighbors

   !> Routine to write data from the array to the buffer
   subroutine array2buffer(ndims, dims, begin, end, array, array_size, buffer, buffer_size, buffer_start_index)
      integer, intent(in) :: ndims, buffer_start_index, array_size, buffer_size
      integer, dimension(ndims), intent(in) :: dims, begin, end
      real, dimension(array_size), intent(in) :: array
      real, dimension(buffer_size), intent(inout) :: buffer

      integer :: ii, jj, global_index, buffer_global_index

      buffer_global_index = 0

      if(ndims == 2) then
         do ii = begin(1), end(1)
            do jj = begin(2), end(2)
               global_index = IDX_XD(ndims, dims, [ii, jj])
               buffer(buffer_start_index + buffer_global_index) = array(global_index)
               buffer_global_index = buffer_global_index + 1
            end do
         end do
      end if

   end subroutine array2buffer

   !> Routine to write data from the buffer to the array
   subroutine buffer2array(ndims, dims, begin, end, array, array_size, buffer, buffer_size, buffer_start_index)
      integer, intent(in) :: ndims, buffer_start_index, array_size, buffer_size
      integer, dimension(ndims), intent(in) :: dims, begin, end
      real, dimension(array_size), intent(inout) :: array
      real, dimension(buffer_size), intent(in) :: buffer

      integer :: ii, jj, global_index, buffer_global_index

      buffer_global_index = 0

      if(ndims == 2) then
         do ii = begin(1), end(1)
            do jj = begin(2), end(2)
               global_index = IDX_XD(ndims, dims, [ii, jj])
               array(global_index) = buffer(buffer_start_index + buffer_global_index)
               buffer_global_index = buffer_global_index + 1
            end do
         end do
      end if

   end subroutine buffer2array

end module rank_module
