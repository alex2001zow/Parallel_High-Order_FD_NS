program main
    use mpi
    use rank_parameters_module
    use utility_functions_module
    implicit none

    type(rank_struct) :: rank_params

    integer :: rank, world_size, ierr

    real, allocatable :: local_matrix(:)

    ! Initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, ierr)
    
    ! Allocate and setup rank parameters
    call allocate_rank_parameters(rank_params, rank, world_size)
  
    ! Print out rank parameters from rank 0
    if(rank_params%rank == 0) THEN
      call print_cartesian_grid(rank_params%ndims, rank_params%processor_dim)
      call print_rank_parameters(rank_params)
    END IF

    allocate(local_matrix(rank_params%total_num_elements))

    deallocate(local_matrix)

    ! Deallocate rank parameters
    call deallocate_rank_parameters(rank_params)
    
    ! Finalize MPI
    call MPI_FINALIZE(ierr)
  end program main
  