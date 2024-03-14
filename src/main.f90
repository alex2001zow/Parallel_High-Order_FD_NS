program main

   call run_simulation()

end program main

subroutine run_simulation()
   use mpi, only: MPI_WTIME, MPI_DOUBLE_PRECISION, MPI_SUM
   use omp_lib
   use constants_module, only: MASTER_RANK, filename_txt, filename_dat
   use utility_functions_module, only: read_input_from_command_line
   use mpi_wrapper_module, only: initialize_mpi_wrapper, finalize_mpi_wrapper, check_openmpi_version, all_reduce_mpi_wrapper
   use rank_module, only: rank_type, create_rank_type, deallocate_rank_type, print_rank_type, write_rank_type_blocks_to_file
   use comm_module, only: print_cartesian_grid
   use solver_module, only: run_solver
   use initialization_module, only: write_initial_condition_and_boundary
   use functions_module, only: Poisson
   implicit none

   integer :: ndims, num_physical_cores
   integer, dimension(:), allocatable :: grid_size, processor_dim

   integer :: rank, world_size
   type(rank_type) :: rank_params

   real, dimension(8) :: result_array_with_timings

   ! Set the number of threads to the number of physical cores. Hopefully not hyperthreading
   num_physical_cores = omp_get_num_procs()
   call omp_set_num_threads(num_physical_cores)

   ! Read input from command line
   call read_input_from_command_line(ndims, grid_size, processor_dim)

   ! Check OpenMPI version
   !call check_openmpi_version()

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size)

   ! Setup rank parameters
   call create_rank_type(ndims, grid_size, processor_dim, int(Poisson,kind=8), rank, world_size, rank_params)

   ! Initialize the block
   call write_initial_condition_and_boundary(rank_params%ndims, rank_params%domain_begin, rank_params%domain_end, &
      rank_params%grid_size, rank_params%block%begin, rank_params%block%size, rank_params%block%matrix, &
      rank_params%FDstencil%dx, rank_params%funcs%initial_condition_func, rank_params%funcs%boundary_condition_func)

   ! Time the program
   result_array_with_timings(5) = MPI_WTIME()

   ! Run the solver
   result_array_with_timings(1:4) = run_solver(rank_params)

   result_array_with_timings(6) = MPI_WTIME()

   result_array_with_timings(7) = result_array_with_timings(6) - result_array_with_timings(5)

   call all_reduce_mpi_wrapper(result_array_with_timings(7), result_array_with_timings(8), 1, &
      int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), rank_params%comm%comm)

   ! Write out the result and timings from the master rank
   if(rank_params%rank == MASTER_RANK) then

      write(*,"(A, E10.3)") "Glob_norm: ", result_array_with_timings(1)
      write(*,"(A, E10.3)") "Rel_norm: ", result_array_with_timings(2)
      write(*,"(A, F10.1)") "Converged: ", result_array_with_timings(3)
      write(*,"(A, F10.1)") "Iterations: ", result_array_with_timings(4)
      write(*,"(A, E10.3, A)") "Total wall time / processors: ", result_array_with_timings(8)/world_size, " seconds"

      ! Write out the cartesian grid from the master rank
      !call print_cartesian_grid(rank_params%comm%comm, rank_params%world_size, rank_params%ndims, filename_txt)
   end if

   ! Write out the rank parameters from each rank. Just for debugging purposes
   call print_rank_type(rank_params, filename_txt)

   ! Write out system solution to a file
   call write_rank_type_blocks_to_file(rank_params, filename_dat)

   ! Deallocate rank parameters
   call deallocate_rank_type(rank_params)

   ! Finalize MPI
   call finalize_mpi_wrapper()

   ! Deallocate grid_size and processor_dim
   if (allocated(grid_size)) deallocate(grid_size)
   if (allocated(processor_dim)) deallocate(processor_dim)

end subroutine run_simulation
