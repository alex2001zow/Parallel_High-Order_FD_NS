module rank_parameters_module
   use mpi_wrapper_module, only: create_cart_communicator_mpi_wrapper, get_cart_coords_mpi_wrapper, &
      cart_rank_mpi_wrapper, change_MPI_COMM_errhandler_mpi_wrapper, &
      original_MPI_COMM_errhandler_mpi_wrapper, isendrecv_mpi_wrapper
   use constants_module, only: neighbor_current_rank, neighbor_non_existant_rank
   use neighbor_types_module, only: get_neighbor_range
   use utility_functions_module, only : IDX_XD
   use initilization_module, only: test_block_structure_2D
   implicit none

   private
   !character(len=MPI_MAX_ERROR_STRING) :: error_string
   !integer :: ierr_string, error_len

   integer :: ierr = 0

   !> Structure to hold the parameters for each rank.
   type rank_struct
      integer :: rank, world_size, cart_comm
      integer :: ndims, num_block_elements, num_block_neighbors, num_sendrecv_elements

      integer, allocatable :: grid_size(:), processor_dim(:), coords(:), neighbors(:), local_size(:)
      integer, allocatable :: neighbor_sendrecv_start_index(:), neighbor_send_request(:), neighbor_recv_request(:)
      real, allocatable :: local_block_matrix(:), neighbor_elements_send(:), neighbor_elements_recv(:)
   end type rank_struct

   public :: rank_struct, setup_rank_parameters, deallocate_rank_parameters, print_rank_parameters
   public :: sendrecv_elements_from_neighbors, buffer2array, array2buffer

contains

   !> Allocate and setup each ranks parameters
   subroutine setup_rank_parameters(rank_val, world_size_val, parameters)
      integer, intent(in) :: rank_val, world_size_val
      type(rank_struct), intent(out) :: parameters

      character(255) :: dim_input_arg, temp_arg

      integer :: ii

      call get_command_argument(1, dim_input_arg)
      read(dim_input_arg,*) parameters%ndims

      parameters%rank = rank_val
      parameters%world_size = world_size_val

      parameters%num_block_neighbors = 3**parameters%ndims

      call allocate_rank_parameters(parameters)

      parameters%num_block_elements = 1
      do ii = 1, parameters%ndims
         call get_command_argument(ii+1, temp_arg)
         read(temp_arg,*) parameters%grid_size(ii)

         call get_command_argument(ii+1+parameters%ndims, temp_arg)
         read(temp_arg,*) parameters%processor_dim(ii)

         ! We have make sure that this is properly divisible
         if(mod(parameters%grid_size(ii), parameters%processor_dim(ii)) /= 0) then
            print *, "Grid size is not divisible by the number of processors in dimension ", ii
            stop
         end if
         parameters%local_size(ii) = parameters%grid_size(ii) / parameters%processor_dim(ii)

         ! The total number of elements we have to allocate per node
         parameters%num_block_elements = parameters%num_block_elements * parameters%local_size(ii)

      end do

      call allocate_block_matrix(parameters)

      ! Create a Cartesian communicator
      call create_cart_communicator_mpi_wrapper(parameters%ndims, parameters%processor_dim, parameters%cart_comm, ierr)

      ! Get the Cartesian coordinates of the current process
      call get_cart_coords_mpi_wrapper(parameters%cart_comm, parameters%rank, parameters%ndims, parameters%coords, ierr)

      call determine_moore_neighbors(parameters)

      call allocate_neighbor_sendrecv_array(parameters)

   end subroutine setup_rank_parameters

   !> Allocate rank parameters. IMPORTANT neighbors are not the correct size yet!
   subroutine allocate_rank_parameters(parameters)
      type(rank_struct), intent(inout) :: parameters

      allocate(parameters%grid_size(parameters%ndims))
      allocate(parameters%processor_dim(parameters%ndims))
      allocate(parameters%coords(parameters%ndims))

      allocate(parameters%neighbors(parameters%num_block_neighbors))
      allocate(parameters%neighbor_sendrecv_start_index(parameters%num_block_neighbors))
      allocate(parameters%neighbor_send_request(parameters%num_block_neighbors))
      allocate(parameters%neighbor_recv_request(parameters%num_block_neighbors))

      allocate(parameters%local_size(parameters%ndims))
   end subroutine allocate_rank_parameters

   !> Deallocate rank parameters
   subroutine deallocate_rank_parameters(parameters)
      type(rank_struct), intent(inout) :: parameters

      ! Deallocate the pointers
      deallocate(parameters%grid_size)
      deallocate(parameters%processor_dim)
      deallocate(parameters%coords)
      deallocate(parameters%neighbors)
      deallocate(parameters%neighbor_sendrecv_start_index)
      deallocate(parameters%local_size)

      call deallocate_block_matrix(parameters)

      deallocate(parameters%neighbor_elements_send)
      deallocate(parameters%neighbor_elements_recv)

      call MPI_Comm_free(parameters%cart_comm, ierr)

   end subroutine deallocate_rank_parameters

   !> Allocate the block matrix
   subroutine allocate_block_matrix(parameters)
      type(rank_struct), intent(inout) :: parameters

      allocate(parameters%local_block_matrix(parameters%num_block_elements))
   end subroutine allocate_block_matrix

   !> Deallocate the block matrix
   subroutine deallocate_block_matrix(parameters)
      type(rank_struct), intent(inout) :: parameters

      deallocate(parameters%local_block_matrix)
   end subroutine deallocate_block_matrix

   !> Print rank parameters
   subroutine print_rank_parameters(parameters)
      type(rank_struct), intent(in) :: parameters
      print *, "Dimensions:", parameters%ndims
      print *, "Rank:", parameters%rank
      print *, "World Size:", parameters%world_size
      print *, "Grid Size:", parameters%grid_size
      print *, "Num Processors per dim:", parameters%processor_dim
      print *, "Cartesian coordinates:", parameters%coords
      print *, "Neighbors:", parameters%neighbors
      print *, "Neighbor sendrecv start index:", parameters%neighbor_sendrecv_start_index
      print *, "Total number of elements to be sent and recieved:", parameters%num_sendrecv_elements
      print *, "Local matrix size:", parameters%local_size
      print *, "Total number of elements in the block:", parameters%num_block_elements
      !print *, "Local block matrix:", parameters%local_block_matrix
   end subroutine print_rank_parameters

   !> Find the neighbors of each rank in a 2D or 3D grid. Would like to make it N-D but that is a bit more complicated.
   subroutine determine_moore_neighbors(parameters)
      type(rank_struct), intent(inout) :: parameters

      integer :: ii, jj, kk, global_index
      integer :: indices(parameters%ndims)
      integer :: rank_of_coords, error

      ! Prevent MPI from aborting when we try to find neighbors that do not exist
      call change_MPI_COMM_errhandler_mpi_wrapper(parameters%cart_comm, ierr)

      ! Determine the ranks of neighboring processes including corners along each dimension. Works for 2D and 3D
      global_index = 1

      if(parameters%ndims == 2) then
         do ii = -1,1,1
            do jj = -1,1,1
               indices = parameters%coords + [ii,jj]

               call cart_rank_mpi_wrapper(parameters%cart_comm, parameters%ndims, indices, rank_of_coords, error, ierr)

               if(error == 1) then
                  parameters%neighbors(global_index) = neighbor_non_existant_rank ! This is a non-existent rank
               else if(ii == 0 .AND. jj == 0) then
                  parameters%neighbors(global_index) = neighbor_current_rank ! This is the current rank
               else
                  parameters%neighbors(global_index) = rank_of_coords
               endif
               global_index = global_index + 1
            end do
         end do
      end if

      if(parameters%ndims == 3) then
         do ii = -1,1,1
            do jj = -1,1,1
               do kk = -1,1,1
                  indices = parameters%coords + [ii,jj,kk]

                  call cart_rank_mpi_wrapper(parameters%cart_comm, parameters%ndims, indices, rank_of_coords, error, ierr)

                  if(error == 1) then
                     parameters%neighbors(global_index) = neighbor_non_existant_rank ! This is a non-existent rank
                  else if(ii == 0 .AND. jj == 0 .AND. kk == 0) then
                     parameters%neighbors(global_index) = neighbor_current_rank ! This is the current rank
                  else
                     parameters%neighbors(global_index) = rank_of_coords
                  endif
                  global_index = global_index + 1
               end do
            end do
         end do
      end if

      ! Restore the original error handler to abort on errors
      call original_MPI_COMM_errhandler_mpi_wrapper(parameters%cart_comm, ierr)

   end subroutine determine_moore_neighbors

   !> Routine to find the elements to be sent and recieved from each neighbor
   subroutine allocate_neighbor_sendrecv_array(parameters)
      type(rank_struct), intent(inout) :: parameters

      integer :: begin(parameters%ndims), end(parameters%ndims)
      integer :: neighbor_index, neighbor_elements_size

      neighbor_elements_size = 1
      do neighbor_index = 1, parameters%num_block_neighbors
         if(parameters%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            call get_neighbor_range(parameters%ndims, neighbor_index, parameters%local_size, begin, end)
            neighbor_elements_size = product(end - begin + 1)
            parameters%neighbor_sendrecv_start_index(neighbor_index) = parameters%num_sendrecv_elements
            parameters%num_sendrecv_elements = parameters%num_sendrecv_elements + neighbor_elements_size
         else
            parameters%neighbor_sendrecv_start_index(neighbor_index) = -1
         endif
      end do

      allocate(parameters%neighbor_elements_send(parameters%num_sendrecv_elements))
      allocate(parameters%neighbor_elements_recv(parameters%num_sendrecv_elements))

   end subroutine allocate_neighbor_sendrecv_array

   !> Routine to send and recieve elements from each neighbor
   subroutine sendrecv_elements_from_neighbors(parameters)
      type(rank_struct), intent(inout) :: parameters

      integer :: begin(parameters%ndims), end(parameters%ndims)
      integer :: neighbor_index, neighbor_elements_size

      do neighbor_index = 1, parameters%num_block_neighbors
         if(parameters%neighbors(neighbor_index) > neighbor_non_existant_rank) then
            call get_neighbor_range(parameters%ndims, neighbor_index, parameters%local_size, begin, end)
            neighbor_elements_size = product(end - begin + 1)
            call isendrecv_mpi_wrapper(parameters%num_sendrecv_elements, parameters%neighbor_elements_send, &
               parameters%neighbor_elements_recv, parameters%neighbor_sendrecv_start_index(neighbor_index), &
               neighbor_elements_size, parameters%neighbors(neighbor_index),&
               parameters%cart_comm, parameters%neighbor_send_request(neighbor_index), &
               parameters%neighbor_recv_request(neighbor_index), ierr)
         end if
      end do

   end subroutine sendrecv_elements_from_neighbors

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
               global_index = IDX_XD(ndims,dims,[jj, ii])
               array(global_index) = buffer(buffer_start_index + buffer_global_index)
               buffer_global_index = buffer_global_index + 1
            end do
         end do
      endif

      if(ndims == 3) then
         do ii = begin(1),end(1)
            do jj = begin(2),end(2)
               do kk = begin(3),end(3)
                  global_index = IDX_XD(ndims,dims,[kk, jj, ii])
                  array(global_index) = buffer(buffer_start_index + buffer_global_index)
                  buffer_global_index = buffer_global_index + 1
               end do
            end do
         end do
      endif

   end subroutine buffer2array

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
               global_index = IDX_XD(ndims, dims, [jj, ii])
               buffer(buffer_start_index + buffer_global_index) = array(global_index)
               buffer_global_index = buffer_global_index + 1
            end do
         end do
      endif

      if(ndims == 3) then
         do ii = begin(1),end(1)
            do jj = begin(2),end(2)
               do kk = begin(3),end(3)
                  global_index = IDX_XD(ndims, dims, [kk, jj, ii])
                  buffer(buffer_start_index + buffer_global_index) = array(global_index)
                  buffer_global_index = buffer_global_index + 1
               end do
            end do
         end do
      endif

   end subroutine array2buffer

end module rank_parameters_module
