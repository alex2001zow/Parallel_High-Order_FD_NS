program main

   call run_simulation()

end program main

subroutine run_simulation()
   use constants_module, only: MASTER_RANK, filename
   use mpi_wrapper_module, only: initialize_mpi_wrapper, finalize_mpi_wrapper
   use rank_module, only: rank_type, create_rank_type, deallocate_rank_type, print_rank_type
   use comm_module, only: print_cartesian_grid
   use solver_module, only: run_solver
   use initialization_module, only: initialize_block_2D
   implicit none

   type(rank_type) :: rank_params
   integer :: rank, world_size
   logical :: converged

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size)

   ! Setup rank parameters
   call create_rank_type(rank, world_size, rank_params)

   ! Initialize the block
   call initialize_block_2D(rank_params%ndims, rank_params%grid_size, rank_params%block%begin,&
      rank_params%block%size, rank_params%block%matrix, rank)

   ! Run the solver
   converged = run_solver(rank_params)

   ! Write out the cartesian grid from the master rank
   if(rank_params%rank == MASTER_RANK) then
      call print_cartesian_grid(rank_params%comm%comm, rank_params%world_size, rank_params%ndims, filename)
   end if

   ! Write out the rank parameters from each rank. Just for debugging purposes
   call print_rank_type(rank_params, filename)

   ! Deallocate rank parameters
   call deallocate_rank_type(rank_params)

   ! Finalize MPI
   call finalize_mpi_wrapper()

end subroutine run_simulation


