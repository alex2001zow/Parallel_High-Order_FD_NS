module rank_module
   use constants_module, only: neighbor_current_rank, neighbor_non_existant_rank
   use functions_module, only: FunctionPair, set_function_pointers
   use mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_PRECISION, MPI_INT
   use comm_module, only : comm_type, create_cart_comm_type, deallocate_cart_comm_type, print_cart_comm_type, &
      get_neighbor_indices
   use block_module, only : block_type, create_block_type, deallocate_block_type, print_block_type
   use FD_module, only : FDstencil_type, create_finite_difference_stencil, &
      deallocate_finite_difference_stencil, print_finite_difference_stencil
   use mpi_wrapper_module, only: create_cart_communicator_mpi_wrapper, get_cart_coords_mpi_wrapper, &
      cart_rank_mpi_wrapper, change_MPI_COMM_errhandler_mpi_wrapper, &
      original_MPI_COMM_errhandler_mpi_wrapper, isendrecv_mpi_wrapper, waitall_mpi_wrapper, free_communicator_mpi_wrapper, &
      write_to_file_mpi_wrapper
   use utility_functions_module, only : IDX_XD, IDX_XD_INV, sleeper_function
   implicit none

   private

   !> Structure to hold the parameters for each rank.
   type rank_type
      integer :: ndims, rank, world_size
      integer, dimension(:), allocatable :: grid_size, processor_dim
      real, dimension(:), allocatable :: domain_begin, domain_end

   end type rank_type

   public :: rank_type
   public :: create_rank_type, deallocate_rank_type
   public :: communicate_step
   public :: print_rank_type, write_rank_type_blocks_to_file

contains

   !> Routine to create the rank_type
   subroutine create_rank_type(ndims, rank, world_size, grid_size, processor_dim, domain_begin, domain_end, parameters)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(ndims), intent(in) :: grid_size, processor_dim
      real, dimension(ndims), intent(in) :: domain_begin, domain_end
      type(rank_type), intent(inout) :: parameters

      parameters%ndims = ndims
      parameters%rank = rank
      parameters%world_size = world_size

      call allocate_rank_type(parameters)
      parameters%grid_size = grid_size
      parameters%processor_dim = processor_dim

      parameters%domain_begin = domain_begin
      parameters%domain_end = domain_end

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

   end subroutine deallocate_rank_type

   !> Routine to print the rank_type
   subroutine print_rank_type(parameters_in, iounit)
      type(rank_type), intent(in) :: parameters_in
      integer, intent(in) :: iounit

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "rank_type: "
      write(iounit, *)

      write(iounit, *) "ndims:", parameters_in%ndims
      write(iounit, *) "rank: ", parameters_in%rank
      write(iounit, *) "world_size: ", parameters_in%world_size

      write(iounit, *) "grid_size: ", parameters_in%grid_size
      write(iounit, *) "processor_dim: ", parameters_in%processor_dim

   end subroutine print_rank_type

   !> Routine to write the block data to a file
   subroutine write_rank_type_blocks_to_file(parameters_in, comm_in, block_in, filename_system_solution_in)
      character(255), intent(in) :: filename_system_solution_in
      type(rank_type), intent(in) :: parameters_in
      type(comm_type), intent(in) :: comm_in
      type(block_type), intent(in) :: block_in

      integer, dimension(parameters_in%ndims) :: block_dims, block_begin
      integer :: offset_solution, elements_to_write, starting_offset, extent
      integer :: num_blocks_type_vector, whole_block_count

      integer :: size_real

      size_real = sizeof(0.0d0) ! Size of a double precision number in bytes

      block_dims = block_in%size ! Size of the block (including ghost points)
      block_begin = block_in%begin ! Start index of the block (including ghost points)

      num_blocks_type_vector = 1
      whole_block_count = product(block_dims)

      elements_to_write = whole_block_count ! Number of elements to write to the file
      offset_solution = (parameters_in%rank) * elements_to_write * size_real ! Where in the file to start writing the true block. Times 8 for double precision
      starting_offset = product(block_begin) * size_real ! The starting offset for the true block in the full block. This is for MPI_TYPE_VECTOR.
      extent = elements_to_write * size_real ! Where the true block ends in the full block. Times 8 for double precision

      call write_to_file_mpi_wrapper(parameters_in%ndims, & ! number of dimensions
         filename_system_solution_in, & ! name of the file
         comm_in%comm, & ! MPI communicator
         num_blocks_type_vector, & ! num blocks for MPI_TYPE_VECTOR. This is the total number of blocks (including ghost points)
         whole_block_count, & ! block_length (no ghost points)
         whole_block_count, & ! full block length (stride) (including ghost points)
         0, & ! starting offset for MPI_TYPE_RESIZED. Currently found using MPI_GET_EXTENT for just the whole block.
         0, & ! extent for the MPI_TYPE_RESIZED. Currently found using MPI_GET_EXTENT for just the whole block.
         parameters_in%grid_size, & ! global dims
         block_dims, & ! local dims of the true block (no ghost points)
         block_begin - 1, & ! start index of the true block (no ghost points). Minus 1 to convert to 0 based indexing for MPI
         block_in%matrix, & ! the block matrix
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

      integer :: global_index,  array_global_index, buffer_global_index
      integer, dimension(ndims) :: local_dims, index

      local_dims = end - begin + 1

      do global_index = 1, product(local_dims)
         index = IDX_XD_INV(ndims, local_dims, global_index)
         index = begin + index - 1

         array_global_index = IDX_XD(ndims, dims, index)
         buffer_global_index = buffer_start_index + global_index-1

         buffer(buffer_global_index) = array(array_global_index)
      end do

   end subroutine array2buffer

   !> Routine to write data from the buffer to the array
   subroutine buffer2array(ndims, dims, begin, end, array, array_size, buffer, buffer_size, buffer_start_index)
      integer, intent(in) :: ndims, buffer_start_index, array_size, buffer_size
      integer, dimension(ndims), intent(in) :: dims, begin, end
      real, dimension(array_size), intent(inout) :: array
      real, dimension(buffer_size), intent(in) :: buffer

      integer :: global_index, array_global_index, buffer_global_index
      integer, dimension(ndims) :: local_dims, index

      local_dims = end - begin + 1

      do global_index = 1, product(local_dims)
         index = IDX_XD_INV(ndims, local_dims, global_index)
         index = begin + index - 1

         array_global_index = IDX_XD(ndims, dims, index)
         buffer_global_index = buffer_start_index + global_index-1

         array(array_global_index) = buffer(buffer_global_index)
      end do

   end subroutine buffer2array

end module rank_module
