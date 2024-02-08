program main
   use constants_module, only: MASTER_RANK, filename
   use mpi_wrapper_module, only: initialize_mpi_wrapper, finalize_mpi_wrapper
   use rank_parameters_module, only: rank_type, setup_rank_parameters, deallocate_rank_parameters
   use solver_module, only: run_solver
   use print_module, only: print_cartesian_grid, print_rank_parameters
   use initialization_module, only: initialize_block_2D
   implicit none

   type(rank_type) :: rank_params
   integer :: rank, world_size
   logical :: converged

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size)

   ! Setup rank parameters
   call setup_rank_parameters(rank, world_size, rank_params)

   ! Initialize the block
   call initialize_block_2D(rank_params%ndims, rank_params%grid_size, rank_params%begin_block,&
      rank_params%block_size, rank_params%block_matrix, rank_params%rank)

   ! Run the solver
   converged = run_solver(rank_params)

   ! Write out the cartesian grid from the master rank
   if(rank_params%rank == MASTER_RANK) THEN
      call print_cartesian_grid(rank_params%cart_comm, rank_params%world_size, rank_params%ndims, filename)
   END IF

   ! Write out the rank parameters from each rank. Just for debugging purposes
   call print_rank_parameters(rank_params, filename)

   ! Deallocate rank parameters
   call deallocate_rank_parameters(rank_params)

   ! Finalize MPI
   call finalize_mpi_wrapper()
end program main

