program main

   call run_simulation()

end program main

subroutine run_simulation()
   use omp_lib
   use constants_module, only: MASTER_RANK, filename_txt, filename_dat
   use utility_functions_module, only: read_input_from_command_line
   use mpi_wrapper_module, only: initialize_mpi_wrapper, finalize_mpi_wrapper, check_openmpi_version
   use rank_module, only: rank_type, create_rank_type, deallocate_rank_type, print_rank_type, write_rank_type_blocks_to_file
   use comm_module, only: print_cartesian_grid
   use solver_module, only: run_solver
   use initialization_module, only: write_initial_condition_and_boundary
   use functions_module, only: Poisson2D
   implicit none

   integer :: ndims, num_physical_cores
   integer, dimension(:), allocatable :: grid_size, processor_dim

   integer :: rank, world_size
   type(rank_type) :: rank_params

   real, dimension(4) :: result_array

   num_physical_cores = omp_get_num_procs()
   call omp_set_num_threads(num_physical_cores)

   ! Read input from command line
   call read_input_from_command_line(ndims, grid_size, processor_dim)

   ! Check OpenMPI version
   !call check_openmpi_version()

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size)

   ! Setup rank parameters
   call create_rank_type(ndims, grid_size, processor_dim, INT(Poisson2D,kind=8), rank, world_size, rank_params)

   ! Initialize the block
   call write_initial_condition_and_boundary(rank_params%ndims, rank_params%domain_begin, &
      rank_params%block%begin, rank_params%block%size, rank_params%block%matrix, &
      rank_params%FDstencil%dx, rank_params%funcs%initial_condition_func, rank_params%funcs%boundary_condition_func)

   ! Run the solver
   result_array = run_solver(rank_params)

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

   ! Deallocate grid_size and processor_dim
   if (allocated(grid_size)) deallocate(grid_size)
   if (allocated(processor_dim)) deallocate(processor_dim)

end subroutine run_simulation
