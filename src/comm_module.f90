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
      integer, allocatable :: coords(:), neighbors(:)
      integer, allocatable :: neighbor_send_request(:), neighbor_recv_request(:)
   end type comm_type

   public :: comm_type, create_cart_comm_type, deallocate_cart_comm_type, print_cart_comm_type

contains

   !> Allocate a comm_type object
   subroutine create_cart_comm_type(ndims, processor_dim, rank, comm_output)
      integer, intent(in) :: ndims, rank
      integer, dimension(ndims), intent(in) :: processor_dim
      type(comm_type), intent(inout) :: comm_output

      comm_output%num_neighbors = 3**ndims

      call allocate_cart_comm_type(ndims, comm_output)

      call create_cart_communicator(ndims, processor_dim, comm_output%comm)

      call get_cart_coords(comm_output%comm, rank, ndims, comm_output%coords)

      call get_cart_neighbors(ndims, comm_output%coords, comm_output%comm, comm_output%neighbors)

      comm_output%neighbor_send_request = MPI_REQUEST_NULL
      comm_output%neighbor_recv_request = MPI_REQUEST_NULL

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

   !> Find the Moore neighborhood of a rank in a cartesian communicator in 2D and 3D. Want to make it N-D but it is complicated.
   subroutine get_cart_neighbors(ndims, coords, cart_comm, neighbors)
      integer, intent(in) :: ndims, cart_comm
      integer, dimension(ndims), intent(in) :: coords

      integer, dimension(3**ndims), intent(inout) :: neighbors

      integer :: ii, jj, kk, global_index
      integer :: indices(ndims)
      integer :: rank_of_coords, error

      ! Prevent MPI from aborting when we try to find neighbors that do not exist
      call change_MPI_COMM_errhandler_mpi_wrapper(cart_comm)

      ! Determine the ranks of neighboring processes including corners along each dimension. Works for 2D and 3D
      global_index = 1

      if(ndims == 2) then
         do ii = -1,1,1
            do jj = -1,1,1
               indices = coords + [jj,ii]

               call cart_rank_mpi_wrapper(cart_comm, ndims, indices, rank_of_coords, error)

               if(error == 1) then
                  neighbors(global_index) = neighbor_non_existant_rank ! This is a non-existent rank
               else if(ii == 0 .AND. jj == 0) then
                  neighbors(global_index) = neighbor_current_rank ! This is the current rank
               else
                  neighbors(global_index) = rank_of_coords
               endif
               global_index = global_index + 1
            end do
         end do
      end if

      if(ndims == 3) then
         do ii = -1,1,1
            do jj = -1,1,1
               do kk = -1,1,1
                  indices = coords + [kk,jj,ii]

                  call cart_rank_mpi_wrapper(cart_comm, ndims, indices, rank_of_coords, error)

                  if(error == 1) then
                     neighbors(global_index) = neighbor_non_existant_rank ! This is a non-existent rank
                  else if(ii == 0 .AND. jj == 0 .AND. kk == 0) then
                     neighbors(global_index) = neighbor_current_rank ! This is the current rank
                  else
                     neighbors(global_index) = rank_of_coords
                  endif
                  global_index = global_index + 1
               end do
            end do
         end do
      end if

      ! Restore the original error handler to abort on errors
      call original_MPI_COMM_errhandler_mpi_wrapper(cart_comm)

   end subroutine get_cart_neighbors

   !> Allocate a comm_type object
   subroutine allocate_cart_comm_type(ndims, comm)
      integer, intent(in) :: ndims
      type(comm_type), intent(inout) :: comm

      allocate(comm%coords(ndims))
      allocate(comm%neighbors(comm%num_neighbors))
      allocate(comm%neighbor_send_request(comm%num_neighbors))
      allocate(comm%neighbor_recv_request(comm%num_neighbors))

   end subroutine allocate_cart_comm_type

   !> Deallocate a comm_type object
   subroutine deallocate_cart_comm_type(comm)
      type(comm_type), intent(inout) :: comm

      call free_communicator_mpi_wrapper(comm%comm)

      if(allocated(comm%coords)) then
         deallocate(comm%coords)
      endif
      if(allocated(comm%neighbors)) then
         deallocate(comm%neighbors)
      endif
      if(allocated(comm%neighbor_send_request)) then
         deallocate(comm%neighbor_send_request)
      endif
      if(allocated(comm%neighbor_recv_request)) then
         deallocate(comm%neighbor_recv_request)
      endif

   end subroutine deallocate_cart_comm_type

   !> Print the contents of a comm_type object
   subroutine print_cart_comm_type(ndims, comm, iounit)
      integer, intent(in) :: ndims
      type(comm_type), intent(in) :: comm
      integer, intent(in) :: iounit

      integer :: ii, jj, global_index

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "comm_type:"
      write(iounit, *)

      write(iounit, *) "comm: ", comm%comm
      write(iounit, *) "coords: ", comm%coords
      write(iounit, *) "neighbor_send_request: ", comm%neighbor_send_request
      write(iounit, *) "neighbor_recv_request: ", comm%neighbor_recv_request

      write(iounit, *) "Neighbor ranks:"
      if(ndims == 2) then
         do ii = 1, 3
            do jj = 1, 3
               global_index = IDX_XD(2, [3,3], [jj, ii])
               write(iounit, '(I5)', advance='no') comm%neighbors(global_index)
            end do
            write(iounit, *)
         end do
      endif

      if(ndims == 3) then
         write(iounit, *) "NOT DONE YET"
      endif

   end subroutine print_cart_comm_type

end module comm_module
