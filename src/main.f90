program main
   use constants_module
   use mpi_wrapper_module, only: initialize_mpi_wrapper, finalize_mpi_wrapper
   use rank_parameters_module, only: rank_struct, setup_rank_parameters, deallocate_rank_parameters
   use solver_module, only: run_solver
   use print_module, only: print_cartesian_grid, print_rank_parameters
   use initialization_module, only: initialize_block_2D
   implicit none

   type(rank_struct) :: rank_params
   integer :: rank, world_size
   logical :: converged
   character(255) :: filename

   ! Set output filename
   filename = "output/output_from_" // repeat(" ", 255 - len_trim("output/output_from_"))

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size)

   ! Setup rank parameters
   call setup_rank_parameters(rank, world_size, rank_params)

   call initialize_block_2D(rank_params%ndims, rank_params%grid_size, rank_params%begin_block,&
      rank_params%block_size, rank_params%block_matrix, rank_params%rank)

   converged = run_solver(rank_params)

   ! Write out the cartesian grid from rank 0
   if(rank_params%rank == MASTER_RANK) THEN
      call print_cartesian_grid(rank_params%cart_comm, rank_params%world_size, rank_params%ndims, filename)
   END IF

   call print_rank_parameters(rank_params, filename)

   ! Deallocate rank parameters
   call deallocate_rank_parameters(rank_params)

   ! Finalize MPI
   call finalize_mpi_wrapper()
end program main

