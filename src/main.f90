program main

   call run_simulation()

end program main

subroutine run_simulation()
   use constants_module, only: MASTER_RANK, filename_txt, filename_dat
   use mpi_wrapper_module, only: initialize_mpi_wrapper, finalize_mpi_wrapper, check_openmpi_version
   use rank_module, only: rank_type, create_rank_type, deallocate_rank_type, print_rank_type, write_rank_type_blocks_to_file
   use comm_module, only: print_cartesian_grid
   use solver_module, only: run_solver
   use initialization_module, only: initialize_block_2D
   use finite_difference_module, only: FDstencil_type, create_finite_difference_stencil, deallocate_finite_difference_stencil
   implicit none

   type(rank_type) :: rank_params
   integer :: rank, world_size
   real, dimension(4) :: result_array

   ! Check OpenMPI version
   !call check_openmpi_version()

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size)

   ! Setup rank parameters
   call create_rank_type(rank, world_size, rank_params)

   ! Initialize the block
   call initialize_block_2D(rank_params%ndims, rank_params%grid_size, rank_params%domain_begin, &
      rank_params%block%begin, rank_params%block%size, rank_params%block%matrix, rank_params%FDstencil%dx, rank)

   ! Run the solver
   !result_array = run_solver(rank_params)

   ! Write out the cartesian grid from the master rank
   if(rank_params%rank == MASTER_RANK) then
      write(*,*) "Global_norm: ", result_array(1), "Relative norm: ", result_array(2), &
         "Converged: ", result_array(3), "Iterations: ", result_array(4)
      call print_cartesian_grid(rank_params%comm%comm, rank_params%world_size, rank_params%ndims, filename_txt)
   end if

   ! Write out the rank parameters from each rank. Just for debugging purposes
   call print_rank_type(rank_params, filename_txt)

   ! Write out the rank blocks to a binary file
   call write_rank_type_blocks_to_file(rank_params, filename_dat)

   ! Deallocate rank parameters
   call deallocate_rank_type(rank_params)

   ! Finalize MPI
   call finalize_mpi_wrapper()

end subroutine run_simulation


