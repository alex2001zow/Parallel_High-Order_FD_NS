module rank_module
   use mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_PRECISION
   use comm_module, only : comm_type, create_cart_comm_type, deallocate_cart_comm_type, print_cart_comm_type
   use block_module, only : block_type, create_block_type, deallocate_block_type, print_block_type
   use finite_difference_module, only : FDstencil_type, create_finite_difference_stencil, &
      deallocate_finite_difference_stencil, print_finite_difference_stencil
   use mpi_wrapper_module, only: create_cart_communicator_mpi_wrapper, get_cart_coords_mpi_wrapper, &
      cart_rank_mpi_wrapper, change_MPI_COMM_errhandler_mpi_wrapper, &
      original_MPI_COMM_errhandler_mpi_wrapper, isendrecv_mpi_wrapper, waitall_mpi_wrapper, free_communicator_mpi_wrapper, &
      write_to_file_mpi_wrapper
   use constants_module, only: neighbor_current_rank, neighbor_non_existant_rank
   use neighbor_types_module, only: get_neighbor_range
   use utility_functions_module, only : IDX_XD, sleeper_function
   implicit none

   private

   !> Structure to hold the parameters for each rank.
   type rank_type
      integer :: ndims, rank, world_size
      integer, dimension(:), allocatable :: grid_size, processor_dim
      real, dimension(:), allocatable :: domain_begin, domain_end

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

      call create_cart_comm_type(parameters%ndims, parameters%processor_dim, parameters%rank, parameters%comm)

      call create_finite_difference_stencil(parameters%ndims, num_derivatives, derivatives, derivatives_order, derivatives_sign,&
         dx, alphas, betas, parameters%FDstencil)

      call create_block_type(parameters%ndims, parameters%comm, parameters%grid_size/parameters%processor_dim,&
         parameters%comm%coords, parameters%block)

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
      endif

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
   subroutine write_rank_type_blocks_to_file(filename, parameters)
      character(255), intent(in) :: filename
      type(rank_type), intent(in) :: parameters

      integer :: start_global_index, offset

      start_global_index = IDX_XD(parameters%ndims, parameters%grid_size, parameters%block%begin)

      offset = (start_global_index-1)*8

      call write_to_file_mpi_wrapper(filename, parameters%comm%comm, parameters%block%matrix, &
         parameters%block%num_elements, INT(MPI_DOUBLE_PRECISION,kind=8), offset)

   end subroutine write_rank_type_blocks_to_file

   !> Routine to communicate with the neighbors
   !! We could edit this routine to initate the recieve first then write to the send buffer, initate send and then check if recv then write to buffer then wait for send and done.
   subroutine communicate_step(parameters)
      type(rank_type), intent(inout) :: parameters

      integer, dimension(parameters%ndims) :: begin, end
      integer :: neighbor_index

      do neighbor_index = 1, parameters%comm%num_neighbors
         if(parameters%comm%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            call get_neighbor_range(parameters%ndims, neighbor_index, parameters%block%size, begin, end)
            call array2buffer(parameters%ndims, parameters%block%size, begin, end, &
               parameters%block%matrix, parameters%block%num_elements, parameters%block%elements_send, &
               parameters%block%num_sendrecv_elements, parameters%block%sendrecv_start_index(neighbor_index))
         end if
      end do

      call sendrecv_elements_from_neighbors(parameters)

      call waitall_mpi_wrapper(parameters%comm%num_neighbors, parameters%comm%neighbor_send_request)
      call waitall_mpi_wrapper(parameters%comm%num_neighbors, parameters%comm%neighbor_recv_request)

      do neighbor_index = 1, parameters%comm%num_neighbors
         if(parameters%comm%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            call get_neighbor_range(parameters%ndims, neighbor_index, parameters%block%size, begin, end)
            call buffer2array(parameters%ndims, parameters%block%size, begin, end, &
               parameters%block%matrix, parameters%block%num_elements, parameters%block%elements_recv, &
               parameters%block%num_sendrecv_elements, parameters%block%sendrecv_start_index(neighbor_index))
         end if
      end do

   end subroutine communicate_step

   !> Routine to send and recieve elements from each neighbor
   subroutine sendrecv_elements_from_neighbors(parameters)
      type(rank_type), intent(inout) :: parameters

      integer, dimension(parameters%ndims) :: begin, end
      integer :: neighbor_index, neighbor_elements_size

      do neighbor_index = 1, parameters%comm%num_neighbors
         if(parameters%comm%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            call get_neighbor_range(parameters%ndims, neighbor_index, parameters%block%size, begin, end)
            neighbor_elements_size = product(end - begin + 1)
            call isendrecv_mpi_wrapper(parameters%block%num_sendrecv_elements, parameters%block%elements_send, &
               parameters%block%elements_recv, parameters%block%sendrecv_start_index(neighbor_index), &
               neighbor_elements_size, parameters%comm%neighbors(neighbor_index),&
               parameters%comm%comm, parameters%comm%neighbor_send_request(neighbor_index), &
               parameters%comm%neighbor_recv_request(neighbor_index))
         else
            parameters%comm%neighbor_send_request(neighbor_index) = MPI_REQUEST_NULL
            parameters%comm%neighbor_recv_request(neighbor_index) = MPI_REQUEST_NULL
         end if
      end do

   end subroutine sendrecv_elements_from_neighbors

   !> Routine to write data from the array to the buffer
   subroutine array2buffer(ndims, dims, begin, end, array, array_size, buffer, buffer_size, buffer_start_index)
      integer, intent(in) :: ndims, buffer_start_index, array_size, buffer_size
      integer, dimension(ndims), intent(in) :: dims, begin, end
      real, dimension(array_size), intent(in) :: array
      real, dimension(buffer_size), intent(inout) :: buffer

      integer :: ii, jj, kk, global_index, buffer_global_index

      buffer_global_index = 0
      if(ndims == 2) then
         do ii = begin(1),end(1)
            do jj = begin(2),end(2)
               global_index = IDX_XD(ndims, dims, [ii, jj])
               buffer(buffer_start_index + buffer_global_index) = array(global_index)
               buffer_global_index = buffer_global_index + 1
            end do
         end do
      endif

      if(ndims == 3) then
         do ii = begin(1),end(1)
            do jj = begin(2),end(2)
               do kk = begin(3),end(3)
                  global_index = IDX_XD(ndims, dims, [ii, jj, kk])
                  buffer(buffer_start_index + buffer_global_index) = array(global_index)
                  buffer_global_index = buffer_global_index + 1
               end do
            end do
         end do
      endif

   end subroutine array2buffer

   !> Routine to write data from the buffer to the array
   subroutine buffer2array(ndims, dims, begin, end, array, array_size, buffer, buffer_size, buffer_start_index)
      integer, intent(in) :: ndims, buffer_start_index, array_size, buffer_size
      integer, dimension(ndims), intent(in) :: dims, begin, end
      real, dimension(array_size), intent(inout) :: array
      real, dimension(buffer_size), intent(in) :: buffer

      integer :: ii, jj, kk, global_index, buffer_global_index

      buffer_global_index = 0
      if(ndims == 2) then
         do ii = begin(1),end(1)
            do jj = begin(2),end(2)
               global_index = IDX_XD(ndims,dims,[ii, jj])
               array(global_index) = buffer(buffer_start_index + buffer_global_index)
               buffer_global_index = buffer_global_index + 1
            end do
         end do
      endif

      if(ndims == 3) then
         do ii = begin(1),end(1)
            do jj = begin(2),end(2)
               do kk = begin(3),end(3)
                  global_index = IDX_XD(ndims,dims,[ii, jj, kk])
                  array(global_index) = buffer(buffer_start_index + buffer_global_index)
                  buffer_global_index = buffer_global_index + 1
               end do
            end do
         end do
      endif

   end subroutine buffer2array

end module rank_module
