program main

   call run_simulation()

end program main

subroutine run_simulation()
   use mpi, only: MPI_WTIME, MPI_DOUBLE_PRECISION, MPI_SUM
   use omp_lib
   use constants_module, only: MASTER_RANK, filename_txt, filename_dat
   use utility_functions_module, only: read_input_from_command_line, open_txt_file, close_txt_file
   use mpi_wrapper_module, only: initialize_mpi_wrapper, finalize_mpi_wrapper, check_openmpi_version, all_reduce_mpi_wrapper

   use rank_module, only: rank_type, create_rank_type, deallocate_rank_type, &
      print_rank_type, write_rank_type_blocks_to_file
   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type, print_cartesian_grid
   use finite_difference_module, only: FDstencil_type, create_finite_difference_stencil, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil
   use block_module, only: block_type, create_block_type, deallocate_block_type, &
      print_block_type
   use functions_module, only: FunctionPair, set_function_pointers, PoissonModel

   use initialization_module, only: write_initial_condition_and_boundary
   use solver_module, only: run_solver
   implicit none

   integer :: ndims, num_physical_cores, iounit
   integer, dimension(:), allocatable :: grid_size, processor_dim
   real, dimension(:), allocatable :: domain_begin, domain_end, dx

   integer :: rank, world_size
   type(rank_type) :: rank_params
   type(comm_type) :: comm_params
   type(FDstencil_type) :: FDstencil_params
   type(block_type) :: block_params
   type(FunctionPair) :: funcs_params

   real, dimension(8) :: result_array_with_timings

   !! Stencil setup. The following is for the Poisson equation with Dirichlet boundary conditions

   ! ! ! 1D-case:
   ! integer, parameter :: num_derivatives = 1
   ! integer, dimension(num_derivatives), parameter :: derivatives = [3]
   ! integer, dimension(num_derivatives), parameter :: derivatives_order = [2]
   ! real, dimension(num_derivatives), parameter :: derivatives_sign = [1]
   ! integer, dimension(1), parameter :: alphas = [1], betas = [1]

   ! 2D-case:
   integer, parameter :: num_derivatives = 2
   integer, dimension(2*num_derivatives), parameter :: derivatives = [1,3,3,1]
   integer, dimension(2*num_derivatives), parameter :: derivatives_order = [0,2,2,0]
   real, dimension(num_derivatives), parameter :: derivatives_sign = [1,1]
   integer, dimension(2), parameter :: alphas = [1,1], betas = [1,1]

   ! ! 3D-case:
   ! integer, parameter :: num_derivatives = 3
   ! integer, dimension(3*num_derivatives), parameter :: derivatives = [1,1,3,1,3,1,3,1,1]
   ! integer, dimension(3*num_derivatives), parameter :: derivatives_order = [0,0,2,0,2,0,2,0,0]
   ! real, dimension(num_derivatives), parameter :: derivatives_sign = [1,1,1]
   ! integer, dimension(3), parameter :: alphas = [1,1,1], betas = [1,1,1]

   ! Set the number of threads to the number of physical cores. Hopefully no hyperthreading
   num_physical_cores = omp_get_num_procs()
   call omp_set_num_threads(num_physical_cores)

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size)

   ! MPI-setup. We can also read input from the command line
   ! call read_input_from_command_line(ndims, grid_size, processor_dim)
   ndims = 2
   domain_begin = [0.0,0.0]
   domain_end = [1.0,1.0]
   grid_size = [8,8]
   processor_dim = [1,1]

   !! Create elements needed for the Poisson equation

   ! Setup rank parameters
   call create_rank_type(ndims, rank, world_size, grid_size, processor_dim, domain_begin, domain_end, rank_params)

   call set_function_pointers(int(PoissonModel,kind=8), funcs_params)

   call create_cart_comm_type(rank_params%ndims, rank_params%grid_size, rank_params%processor_dim, rank_params%rank, comm_params)

   dx = abs((rank_params%domain_end - rank_params%domain_begin)) / (rank_params%grid_size - 1)

   call create_finite_difference_stencil(rank_params%ndims, num_derivatives, derivatives, derivatives_order,&
      derivatives_sign, dx, alphas, betas, FDstencil_params)

   call create_block_type(rank_params%ndims, comm_params, block_params)

   ! Initialize the block
   call write_initial_condition_and_boundary(rank_params%ndims, rank_params%domain_begin, rank_params%domain_end, &
      rank_params%grid_size, block_params%begin, block_params%size, block_params%matrix, &
      FDstencil_params%dx, funcs_params%initial_condition_func, funcs_params%boundary_condition_func)

   ! Time the program
   result_array_with_timings(5) = MPI_WTIME()

   ! Run the solver
   result_array_with_timings(1:4) = run_solver(rank_params, comm_params, block_params, FDstencil_params, funcs_params)

   result_array_with_timings(6) = MPI_WTIME()

   result_array_with_timings(7) = result_array_with_timings(6) - result_array_with_timings(5)

   call all_reduce_mpi_wrapper(result_array_with_timings(7), result_array_with_timings(8), 1, &
      int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_params%comm)

   ! Write out the result and timings from the master rank
   if(rank_params%rank == MASTER_RANK) then

      write(*,"(A, E10.3)") "Glob_norm: ", result_array_with_timings(1)
      write(*,"(A, E10.3)") "Rel_norm: ", result_array_with_timings(2)
      write(*,"(A, F10.1)") "Converged: ", result_array_with_timings(3)
      write(*,"(A, F10.1)") "Iterations: ", result_array_with_timings(4)
      write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(8)/world_size, " seconds"

      ! Write out the cartesian grid from the master rank
      call print_cartesian_grid(rank_params%ndims, comm_params%comm, rank_params%world_size, rank_params%processor_dim,&
         filename_txt)
   end if

   ! Write out our setup to a file. Just for debugging purposes
   call open_txt_file(filename_txt, rank_params%rank, iounit)

   call print_rank_type(rank_params, iounit)

   call print_cart_comm_type(rank_params%ndims, comm_params, iounit)

   call print_finite_difference_stencil(rank_params%ndims, FDstencil_params, iounit)

   call print_block_type(rank_params%ndims, block_params, iounit)

   call close_txt_file(iounit)

   ! Write out system solution to a file
   call write_rank_type_blocks_to_file(rank_params, comm_params, block_params, filename_dat)

   ! Deallocate data
   call deallocate_rank_type(rank_params)

   call deallocate_cart_comm_type(comm_params)

   call deallocate_finite_difference_stencil(FDstencil_params)

   call deallocate_block_type(block_params)

   ! Finalize MPI
   call finalize_mpi_wrapper()

   ! Deallocate the system setup
   if (allocated(grid_size)) deallocate(grid_size)
   if (allocated(processor_dim)) deallocate(processor_dim)
   if (allocated(domain_begin)) deallocate(domain_begin)
   if (allocated(domain_end)) deallocate(domain_end)
   if (allocated(dx)) deallocate(dx)

end subroutine run_simulation
