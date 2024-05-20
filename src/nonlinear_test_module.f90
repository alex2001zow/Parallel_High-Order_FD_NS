module nonlinear_test_module
   use constants_module, only : pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type, print_cartesian_grid
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, apply_FDstencil
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type
   use functions_module, only: FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, sleeper_function
   use initialization_module, only: write_initial_condition_and_boundary

   use solver_module, only: SolverPtrType, set_SystemSolver_pointer, SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, choose_iterative_solver
   implicit none

   private

   public :: nonlinear_1D_test_main

contains

   subroutine nonlinear_1D_test_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      integer :: ndims
      integer, dimension(:), allocatable :: grid_size, processor_dims, stencil_sizes
      real, dimension(:), allocatable :: domain_begin, domain_end
      type(SolverParamsType) :: solver_params

      ndims = 1

      allocate(grid_size(ndims))
      allocate(processor_dims(ndims))
      allocate(domain_begin(ndims))
      allocate(domain_end(ndims))
      allocate(stencil_sizes(ndims))

      grid_size = 8
      processor_dims = 2
      domain_begin = 0
      domain_end = 1
      stencil_sizes = 3

      ! Set the solver parameters tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(1e-6, 1e-2, 10000, 2, solver_params)

      call nonlinear_1D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         stencil_sizes, solver_params)

      ! Deallocate the system setup
      if (allocated(grid_size)) deallocate(grid_size)
      if (allocated(processor_dims)) deallocate(processor_dims)
      if (allocated(domain_begin)) deallocate(domain_begin)
      if (allocated(domain_end)) deallocate(domain_end)

   end subroutine nonlinear_1D_test_main

   subroutine nonlinear_1D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      stencil_sizes, solver_params)
      integer, intent(in) :: rank, world_size
      integer, dimension(1), intent(in) :: grid_size, processor_dims, stencil_sizes
      real, dimension(1), intent(in) :: domain_begin, domain_end
      type(SolverParamsType), intent(in) :: solver_params

      integer, parameter :: num_derivatives = 2
      integer, dimension(num_derivatives), parameter :: derivatives = [1,2]

      call nonlinear_1D_test(1, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, stencil_sizes, solver_params)

   end subroutine nonlinear_1D_analytical

   subroutine nonlinear_1D_test(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      num_derivatives, derivatives, stencil_sizes, solver_params)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(:), intent(in) :: grid_size, processor_dims
      real, dimension(:), intent(in) :: domain_begin, domain_end

      integer, intent(in) :: num_derivatives
      integer, dimension(:), intent(in) :: derivatives, stencil_sizes
      type(SolverParamsType), intent(in) :: solver_params

      integer, dimension(5) :: num_data_elements

      integer, dimension(ndims) ::  bc_begin, bc_end, ghost_begin, ghost_end, stencil_begin, stencil_end
      integer, dimension(ndims) ::  begin, end
      integer :: iounit
      real, dimension(4) :: result_array_with_timings

      type(comm_type) :: comm_params
      type(FDstencil_type) :: FDstencil_params, FDstencil_params2
      type(block_type) :: block_params
      type(FunctionPair) :: funcs_params
      type(SolverPtrType) :: SystemSolver
      type(ResultType) :: result

      num_data_elements = 1

      bc_begin = 0
      bc_end = 0
      ghost_begin = stencil_sizes/2
      ghost_end = stencil_sizes/2
      stencil_begin = stencil_sizes/2
      stencil_end = stencil_sizes/2

      call set_SystemSolver_pointer(nonlinear_1D_solve_system, SystemSolver)

      call set_function_pointers(num_data_elements, initial_nonlinear_test, boundary_nonlinear_test, &
         f_nonlinear_test, f_nonlinear_test, u_nonlinear_test, funcs_params)

      call create_cart_comm_type(ndims, processor_dims, rank, world_size, comm_params)

      !call sleeper_function(1)

      call create_block_type(ndims, num_data_elements(1), num_data_elements(1), domain_begin, domain_end, grid_size, comm_params, &
         bc_begin, bc_end, ghost_begin, ghost_end, stencil_begin, stencil_end, block_params)

      call create_finite_difference_stencils(block_params%ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      call create_finite_difference_stencils(block_params%ndims, num_derivatives, derivatives, stencil_sizes+2, FDstencil_params2)

      ! Initialize the block
      call write_initial_condition_and_boundary(ndims, num_data_elements(1), domain_begin, domain_end, &
         grid_size, block_params%global_begin_c+1, block_params%local_size, block_params%matrix, &
         block_params%extended_grid_dx, funcs_params)

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      ! Run the solver
      call choose_iterative_solver(comm_params, block_params, FDstencil_params, &
         funcs_params, SystemSolver, solver_params, result)

      if(rank == MASTER_RANK) then
         call print_resultType(result)
      end if

      call choose_iterative_solver(comm_params, block_params, FDstencil_params2, &
         funcs_params, SystemSolver, solver_params, result)

      result_array_with_timings(2) = MPI_WTIME()

      result_array_with_timings(3) = result_array_with_timings(2) - result_array_with_timings(1)

      call all_reduce_mpi_wrapper(result_array_with_timings(3), result_array_with_timings(4), 1, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_params%comm)

      ! Write out the result and timings from the master rank
      if(rank == MASTER_RANK) then
         call print_resultType(result)
         write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(4)/world_size, " seconds"

         ! Write out the cartesian grid from the master rank
         call print_cartesian_grid(ndims, comm_params%comm, world_size, processor_dims, "output/output_from_")
      end if

      ! Write out our setup to a file. Just for debugging purposes
      call open_txt_file("output/output_from_", rank, iounit)

      call print_cart_comm_type(comm_params, iounit)

      call print_finite_difference_stencil(ndims, FDstencil_params, iounit)

      call print_block_type(ndims, block_params, iounit)

      call close_txt_file(iounit)

      ! Write out system solution to a file
      call write_block_data_to_file(block_params%data_layout, "output/system_solution.dat", comm_params%comm, block_params%matrix)

      ! Deallocate data

      call deallocate_cart_comm_type(comm_params)

      call deallocate_finite_difference_stencil(FDstencil_params)

      call deallocate_finite_difference_stencil(FDstencil_params2)

      call deallocate_block_type(block_params)

   end subroutine nonlinear_1D_test

   pure subroutine nonlinear_1D_solve_system(ndims, num_elements, stencil_size, alphas, num_derivatives, &
      stencils, combined_stencils, dims, index, matrix, f_value, F, J)
      integer, intent(in) :: ndims, num_elements, num_derivatives
      integer, dimension(ndims), intent(in) :: stencil_size, alphas, dims, index
      real, dimension(:), target, intent(inout) :: stencils
      real, dimension(:), intent(in) :: matrix
      real, dimension(:), intent(in) :: f_value

      real, dimension(product(stencil_size)), intent(inout) :: combined_stencils
      real, dimension(num_elements), intent(inout) :: F, J

      integer ::num_stencil_elements, center_coefficient_global_index

      real :: u_x, u_xx, u_x_j, u_xx_j

      ! Get the number of stencil elements
      num_stencil_elements = product(stencil_size)

      ! Evaluate the stencils
      call apply_FDstencil(ndims, 1, 0, stencil_size, alphas, stencils(1:num_stencil_elements), &
         dims, index, matrix, u_x)
      call apply_FDstencil(ndims, 1, 0, stencil_size, alphas, stencils(num_stencil_elements+1:2*num_stencil_elements), &
         dims, index, matrix, u_xx)

      F(1) = u_x*u_xx - f_value(1)

      ! Evaluate the Jacobian. This is simply finding the center coefficient
      call IDX_XD(ndims, stencil_size, alphas + 1, center_coefficient_global_index)

      ! Evaluate the Jacobian
      u_x_j = stencils(center_coefficient_global_index)
      u_xx_j = stencils(num_stencil_elements+center_coefficient_global_index)

      J(1) = u_x_j*u_xx + u_x*u_xx_j

   end subroutine nonlinear_1D_solve_system

   pure subroutine initial_nonlinear_test(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      func_val = 20.0
      !call u_nonlinear_test(ndims, global_begin_indices, local_indices, global_domain_begin, global_domain_end, &
      !   global_domain_size, dx, func_val)
      !func_val = func_val + 0.0001

   end subroutine initial_nonlinear_test

   pure subroutine boundary_nonlinear_test(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      integer, dimension(ndims) :: block_index
      logical :: at_boundary

      block_index = global_begin_indices + local_indices - 1

      at_boundary = any(block_index == 1 .or. block_index == global_domain_size)

      if (at_boundary) then
         call u_nonlinear_test(ndims, num_elements, global_begin_indices, local_indices, &
            global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      end if

   end subroutine boundary_nonlinear_test

   pure subroutine f_nonlinear_test(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      real, dimension(ndims) :: point

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      func_val = -pi*pi*pi*product(cos(pi*point)*sin(pi*point))

   end subroutine f_nonlinear_test

   pure subroutine u_nonlinear_test(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      real, dimension(ndims) :: point

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      func_val = product(sin(pi*point))

   end subroutine u_nonlinear_test

end module nonlinear_test_module




