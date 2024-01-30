module rank_parameters_module
  use mpi
  use utility_functions_module
  implicit none

  !> Structure to hold the parameters for each rank. Neighbors are stored in the following order: up down left right front back
  type rank_struct
    integer :: ndims, total_num_elements
    integer :: rank, world_size, cart_comm, ierr
    integer, allocatable :: grid_size(:), processor_dim(:), coords(:), neighbors(:), local_size(:)
    logical, allocatable :: periods(:)
  end type rank_struct

  contains

  !> Allocate and setup each ranks parameters
  subroutine allocate_rank_parameters(parameters, rank_val, world_size_val)
    character(255) :: dim_input_arg, temp_arg
    integer, intent(in) :: rank_val, world_size_val
    type(rank_struct), intent(out) :: parameters

    integer :: ii, global_index

    call get_command_argument(1, dim_input_arg)
    read(dim_input_arg,*) parameters%ndims

    parameters%rank = rank_val

    parameters%world_size = world_size_val

    allocate(parameters%grid_size(parameters%ndims))
    allocate(parameters%processor_dim(parameters%ndims))
    allocate(parameters%periods(parameters%ndims))
    allocate(parameters%coords(parameters%ndims))
    allocate(parameters%neighbors(parameters%ndims*2))
    allocate(parameters%local_size(parameters%ndims))
   
    parameters%total_num_elements = 1
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

      ! We do not have any period in any of the dimensions in the cartesian processor system
      parameters%periods(ii) = .FALSE.

      ! The total number of elements we have to allocate per node
      parameters%total_num_elements = parameters%total_num_elements * parameters%local_size(ii)
      
    end do

    ! Create a Cartesian communicator
    call MPI_Cart_create(MPI_COMM_WORLD, parameters%ndims, parameters%processor_dim, parameters%periods, .TRUE., &
    parameters%cart_comm, parameters%ierr)

    ! Get the Cartesian coordinates of the current process
    call MPI_Cart_coords(parameters%cart_comm, parameters%rank, parameters%ndims, parameters%coords, parameters%ierr)

    ! Determine the ranks of neighboring processes along each dimension
    do ii = 1, parameters%ndims
      global_index = IDX_2D(parameters%ndims, ii, 1)
      call MPI_Cart_shift(parameters%cart_comm, ii-1, 1, parameters%neighbors(global_index),&
      parameters%neighbors(global_index+1), parameters%ierr)
    end do

    ! Determine the corners where neighbors exists
    do ii = 1,  parameters%ndims*2



  end subroutine allocate_rank_parameters

  !> Deallocate rank parameters
  subroutine deallocate_rank_parameters(parameters)
    type(rank_struct), intent(inout) :: parameters
  
    ! Deallocate the pointers
    deallocate(parameters%grid_size)
    deallocate(parameters%processor_dim)
    deallocate(parameters%periods)
    deallocate(parameters%coords)
    deallocate(parameters%neighbors)
    deallocate(parameters%local_size)

    call MPI_Comm_free(parameters%cart_comm, parameters%ierr)
     
  end subroutine deallocate_rank_parameters
  
  !> Print rank parameters
  subroutine print_rank_parameters(parameters)
    type(rank_struct), intent(in) :: parameters
    print *, "Dimensions:", parameters%ndims
    print *, "Rank:", parameters%rank
    print *, "World Size:", parameters%world_size
    print *, "Grid Size:", parameters%grid_size
    print *, "Num Processors per dim:", parameters%processor_dim
    print *, "Periodic:", parameters%periods
    print *, "Cartesian coordinates:", parameters%coords
    print *, "Neighbors:", parameters%neighbors
    print *, "Local matrix size:", parameters%local_size
    print *, "Total number of elements:", parameters%total_num_elements
  end subroutine print_rank_parameters

end module rank_parameters_module
