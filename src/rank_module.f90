module rank_module
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
   use constants_module, only: neighbor_current_rank, neighbor_non_existant_rank
   use utility_functions_module, only : IDX_XD, sleeper_function
   implicit none

   private

   !> Structure to hold the parameters for each rank.
   type rank_type
      integer :: ndims, rank, world_size
      integer, allocatable :: grid_size(:), processor_dim(:)
      real, allocatable :: domain_begin(:), domain_end(:)

      type(comm_type) :: comm
      type(FDstencil_type) :: FDstencil
      type(block_type) :: block

   end type rank_type

   public :: rank_type, create_rank_type, deallocate_rank_type, print_rank_type
   public :: communicate_step
   public :: write_rank_type_blocks_to_file


contains

   !> Routine to create the rank_type
   subroutine create_rank_type(rank, world_size, parameters)
      integer, intent(in) :: rank, world_size
      type(rank_type), intent(out) :: parameters

      character(255) :: dim_input_arg, temp_arg

      integer :: ii

      integer, parameter :: num_derivatives = 2
      integer, dimension(2*num_derivatives), parameter :: derivatives = [1,3,3,1]
      integer, dimension(2*num_derivatives), parameter :: derivatives_order = [0,2,2,0]
      real, dimension(num_derivatives), parameter :: derivatives_sign = [1,1]

      integer, dimension(2), parameter :: alphas = [1,1], betas = [1,1]
      real, dimension(2) :: dx

      call get_command_argument(1, dim_input_arg)
      read(dim_input_arg,*) parameters%ndims

      parameters%rank = rank
      parameters%world_size = world_size

      call allocate_rank_type(parameters)

      do ii = 1, parameters%ndims
         call get_command_argument(ii+1, temp_arg)
         read(temp_arg,*) parameters%grid_size(ii)

         call get_command_argument(ii+1+parameters%ndims, temp_arg)
         read(temp_arg,*) parameters%processor_dim(ii)

         ! We have make sure that the inputs are properly divisible
         if(mod(parameters%grid_size(ii), parameters%processor_dim(ii)) /= 0) then
            print *, "Grid size is not divisible by the number of processors in dimension ", ii
            stop
         end if

         parameters%domain_begin(ii) = 0.0
         parameters%domain_end(ii) = 1.0

         dx(ii) = abs((parameters%domain_end(ii) - parameters%domain_begin(ii))) / (parameters%grid_size(ii)-1)
      end do

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
   subroutine write_rank_type_blocks_to_file(parameters, filename_system_solution, filename_system_solution_layout)
      character(255), intent(in) :: filename_system_solution, filename_system_solution_layout
      type(rank_type), intent(in) :: parameters

      integer :: offset_solution, offset_layout, count_layout
      integer, dimension(1+parameters%ndims*3+1) :: data_layout

      real, dimension(product((parameters%block%size - parameters%comm%offset_end) - &
         (1 - parameters%comm%offset_begin) + 1)) :: true_matrix

      call extract_true_block(parameters%ndims, parameters%block%size, 1 - parameters%comm%offset_begin,&
         parameters%block%size - parameters%comm%offset_end, parameters%block%matrix, true_matrix)

      offset_solution = (parameters%rank) * SIZE(true_matrix) * 8

      ! Layout of the data
      data_layout(:) = [parameters%rank, parameters%comm%coords, parameters%grid_size, parameters%processor_dim, SIZE(true_matrix)]

      count_layout = SIZE(data_layout)

      offset_layout = (parameters%rank) * count_layout * 8

      call write_to_file_mpi_wrapper(filename_system_solution, filename_system_solution_layout,&
         parameters%comm%comm, true_matrix, SIZE(true_matrix), offset_solution, data_layout, count_layout, offset_layout)
   end subroutine write_rank_type_blocks_to_file

   !> Routine to extract the true block from the matrix
   subroutine extract_true_block(ndims, dims, begin_with_offset, end_with_offset, matrix, true_matrix)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims, begin_with_offset, end_with_offset
      real, dimension(product(dims)), intent(in) :: matrix
      real, dimension(product(end_with_offset - begin_with_offset + 1)), intent(out) :: true_matrix

      integer :: ii, jj, matrix_global_index, true_matrix_global_index

      integer, dimension(ndims) :: index

      index(:) = 1

      if(ndims == 2) then
         do ii = begin_with_offset(1), end_with_offset(1)
            do jj = begin_with_offset(2), end_with_offset(2)
               matrix_global_index = IDX_XD(ndims, dims, [ii, jj])
               true_matrix_global_index = IDX_XD(ndims, (end_with_offset - begin_with_offset + 1), index)
               true_matrix(true_matrix_global_index) = matrix(matrix_global_index)
               index(2) = index(2) + 1
            end do
            index(1) = index(1) + 1
            index(2) = 1
         end do
      end if

   end subroutine extract_true_block

   !> Routine to communicate with the neighbors
   !! We could edit this routine to initate the recieve first then write to the send buffer, initate send and then check if recv then write to buffer then wait for send and done.
   subroutine communicate_step(ndims, comm, block)
      integer, intent(in) :: ndims
      type(comm_type), intent(inout) :: comm
      type(block_type), intent(inout) :: block

      integer, dimension(ndims) :: begin, end
      integer :: neighbor_index

      ! Write the data to the send buffer
      do neighbor_index = 1, comm%num_neighbors
         if(comm%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            call get_neighbor_indices(ndims, neighbor_index, comm%neighbor_send_indices, begin, end)
            call array2buffer(ndims, block%size, begin, end, block%matrix, block%num_elements, block%elements_send, &
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
            call buffer2array(ndims, block%size, begin, end, block%matrix, block%num_elements, block%elements_recv, &
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
