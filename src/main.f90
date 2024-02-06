program main
   use constants_module
   use mpi_wrapper_module
   use rank_parameters_module
   use utility_functions_module
   use solver_module
   use print_module
   use initilization_module
   implicit none

   type(rank_struct) :: rank_params
   integer :: rank, world_size, ierr
   logical :: converged
   character(255) :: filename

   ! Set output filename
   filename = "output/output_from_" // repeat(" ", 255 - len_trim("output/output_from_"))

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size, ierr)

   ! Setup rank parameters
   call setup_rank_parameters(rank, world_size, rank_params)

   call initialize_block_2D(rank_params%ndims, rank_params%grid_size, rank_params%begin_block,&
      rank_params%block_size, rank_params%block_matrix, rank_params%rank)

   converged = run_solver(rank_params)

   ! Write out the cartesian grid from rank 0
   if(rank_params%rank == MASTER_RANK) THEN
      call print_cartesian_grid(rank_params%ndims, rank_params%processor_dim, filename)
   END IF

   call print_rank_parameters(rank_params, filename)

   ! Deallocate rank parameters
   call deallocate_rank_parameters(rank_params)

   ! Finalize MPI
   call finalize_mpi_wrapper(ierr)
end program main

