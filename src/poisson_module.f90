module Poisson_module
   use constants_module, only : pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type, print_cartesian_grid
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, apply_FDstencil, update_value_from_stencil, get_FD_coefficients_from_index, &
      calculate_scaled_coefficients
   use block_module, only: block_type, create_block_type, deallocate_block_type, sendrecv_data_neighbors, print_block_type
   use functions_module, only: FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, sleeper_function
   use initialization_module, only: write_initial_condition_and_boundary

   use solver_module, only: SolverPtrType, set_SystemSolver_pointer, SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, choose_iterative_solver, check_convergence
   implicit none

   private

   public :: Poisson_main

contains
   !> The main Poisson subroutine to solve the system
   subroutine Poisson_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      integer :: ndims
      integer, dimension(:), allocatable :: grid_size, processor_dims, stencil_sizes
      real, dimension(:), allocatable :: domain_begin, domain_end
      type(SolverParamsType) :: solver_params

      !call sleeper_function(1)

      ndims = 2

      allocate(grid_size(ndims))
      allocate(processor_dims(ndims))
      allocate(domain_begin(ndims))
      allocate(domain_end(ndims))
      allocate(stencil_sizes(ndims))

      grid_size = 32
      processor_dims = 1
      domain_begin = 0
      domain_end = 1
      stencil_sizes = 3

      ! Set the solver parameters tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(1e-6, 1e-2, 100000, 2, solver_params)

      if(ndims == 1) then
         call Poission_1D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
            stencil_sizes,solver_params)
      end if
      if(ndims == 2) then
         call Poission_2D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
            stencil_sizes, solver_params)
      end if
      if(ndims == 3) then
         call Poission_3D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
            stencil_sizes,solver_params)
      end if

      ! Deallocate the system setup
      if (allocated(grid_size)) deallocate(grid_size)
      if (allocated(processor_dims)) deallocate(processor_dims)
      if (allocated(domain_begin)) deallocate(domain_begin)
      if (allocated(domain_end)) deallocate(domain_end)

   end subroutine Poisson_main

   !> The analytical solution of the 1D Poisson equation
   subroutine Poission_1D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      stencil_sizes, solver_params)
      integer, intent(in) :: rank, world_size
      integer, dimension(1), intent(in) :: grid_size, processor_dims, stencil_sizes
      real, dimension(1), intent(in) :: domain_begin, domain_end
      type(SolverParamsType), intent(in) :: solver_params

      integer, parameter :: num_derivatives = 1
      integer, dimension(num_derivatives), parameter :: derivatives = [2]

      call Poission_analytical_wrapper(1, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, stencil_sizes, solver_params)

   end subroutine Poission_1D_analytical

   !> The analytical solution of the 2D Poisson equation
   subroutine Poission_2D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      stencil_sizes, solver_params)
      integer, intent(in) :: rank, world_size
      integer, dimension(2), intent(in) :: grid_size, processor_dims, stencil_sizes
      real, dimension(2), intent(in) :: domain_begin, domain_end
      type(SolverParamsType), intent(in) :: solver_params

      integer, parameter :: num_derivatives = 2
      integer, dimension(2*num_derivatives), parameter :: derivatives = [0,2,2,0]

      call Poission_analytical_wrapper(2, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, stencil_sizes, solver_params)

   end subroutine Poission_2D_analytical

   !> The analytical solution of the 3D Poisson equation
   subroutine Poission_3D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      stencil_sizes, solver_params)
      integer, intent(in) :: rank, world_size
      integer, dimension(3), intent(in) :: grid_size, processor_dims, stencil_sizes
      real, dimension(3), intent(in) :: domain_begin, domain_end
      type(SolverParamsType), intent(in) :: solver_params

      integer, parameter :: num_derivatives = 3
      integer, dimension(3*num_derivatives), parameter :: derivatives = [0,0,2,0,2,0,2,0,0]

      call Poission_analytical_wrapper(3, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, stencil_sizes, solver_params)

   end subroutine Poission_3D_analytical

   !> The analytical solution of the 2D Poisson equation
   subroutine Poission_analytical_wrapper(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      num_derivatives, derivatives, stencil_sizes, solver_params)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(:), intent(in) :: grid_size, processor_dims
      real, dimension(:), intent(in) :: domain_begin, domain_end

      integer, intent(in) :: num_derivatives
      integer, dimension(:), intent(in) :: derivatives, stencil_sizes
      type(SolverParamsType), intent(in) :: solver_params

      integer, dimension(5) :: num_data_elements

      integer, dimension(ndims) ::  bc_begin, bc_end, ghost_begin, ghost_end, stencil_begin, stencil_end, begin, end
      integer :: iounit
      real, dimension(4) :: result_array_with_timings

      type(comm_type) :: comm_params
      type(FDstencil_type) :: FDstencil_params
      type(block_type) :: block_params
      type(FunctionPair) :: funcs_params
      type(SolverPtrType) :: SystemSolver
      type(ResultType) :: result

      bc_begin = 0
      bc_end = 0
      ghost_begin = stencil_sizes/2
      ghost_end = stencil_sizes/2
      stencil_begin = stencil_sizes/2
      stencil_end = stencil_sizes/2

      num_data_elements = [1,1,1,1,1]

      call set_SystemSolver_pointer(Poisson_solve_system, SystemSolver)

      call set_function_pointers(num_data_elements, initial_poisson, boundary_poisson, &
         f_analytical_poisson, df_analytical_poisson, u_analytical_poisson, funcs_params)

      call create_cart_comm_type(ndims, processor_dims, rank, world_size, comm_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         bc_begin, bc_end, ghost_begin, ghost_end, stencil_begin, stencil_end, block_params)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      ! Initialize the block
      call write_initial_condition_and_boundary(ndims, num_data_elements(1), domain_begin, domain_end, &
         block_params%total_grid_size, block_params%global_begin_c+1, block_params%extended_global_dims-1, &
         block_params%matrix, block_params%dx, funcs_params)

      call sendrecv_data_neighbors(comm_params%comm, block_params, block_params%matrix_ptr)

      ! Set the begin and end for the solver
      begin = 2
      end = block_params%block_dims

      !call sleeper_function(1)

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      ! Run the solver
      !call choose_iterative_solver(comm_params, block_params, FDstencil_params, &
      !   funcs_params, SystemSolver, solver_params, &
      !   begin, end, result)

      ! call GS_Method(comm_params, block_params, FDstencil_params, funcs_params, SystemSolver, solver_params, &
      !    begin, end)

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

      call deallocate_block_type(block_params)

   end subroutine Poission_analytical_wrapper

   !> Solve the Poisson system
   pure subroutine Poisson_solve_system(ndims, num_elements, stencil_size, alphas, num_derivatives, stencils, combined_stencils, &
      dims, index, matrix, f_value, F, J)
      integer, intent(in) :: ndims, num_elements, num_derivatives
      integer, dimension(ndims), intent(in) :: stencil_size, alphas, dims, index
      real, dimension(:), target, intent(inout) :: stencils
      real, dimension(:), intent(in) :: matrix
      real, dimension(:), intent(in) :: f_value

      real, dimension(product(stencil_size)), intent(inout) :: combined_stencils
      real, dimension(num_elements), intent(inout) :: F, J

      integer :: global_index, stencil_begin, stencil_end, num_stencil_elements, center_coefficient_global_index

      combined_stencils = 0.0

      ! Get the number of stencil elements
      num_stencil_elements = product(stencil_size)

      ! Combined stencil into one array since it is linear.
      do global_index = 1, num_derivatives
         stencil_begin = (global_index - 1) * num_stencil_elements + 1
         stencil_end = global_index * num_stencil_elements

         combined_stencils = combined_stencils + stencils(stencil_begin:stencil_end)
      end do

      ! Evaluate the stencil
      call apply_FDstencil(ndims, 1, 0, stencil_size, alphas, combined_stencils, dims, index, matrix, F(1))
      F(1) = F(1) - f_value(1)

      ! Evaluate the Jacobian. This is simply finding the center coefficient
      call IDX_XD(ndims, stencil_size, alphas + 1, center_coefficient_global_index)
      J(1) = combined_stencils(center_coefficient_global_index)

   end subroutine Poisson_solve_system

   subroutine GS_Method(comm_params, block_params, FDstencil_params, funcs_params, SystemSolver, solver_params, &
      begin, end)
      type(comm_type), intent(in) :: comm_params
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), target, intent(inout) :: FDstencil_params
      type(FunctionPair), intent(in) :: funcs_params
      type(SolverPtrType), intent(in) :: SystemSolver
      type(SolverParamsType), intent(in) :: solver_params
      integer, dimension(:), intent(in) :: begin, end

      integer :: ndims, num_elements, num_stencil_elements, it, ii, jj, global_index, converged
      real, dimension(comm_params%ndims) :: global_domain_begin, global_domain_end
      integer, dimension(comm_params%ndims) :: start_dims, dims, index, alphas, global_domain_size, stencil_size
      real, dimension(4) :: norm_array
      real, dimension(:), pointer :: coefficients, dfxx, dfyy
      real :: old_val, new_val
      real, dimension(1) :: f_val
      real, dimension(product(FDstencil_params%stencil_sizes)) :: combined_stencils

      stencil_size = FDstencil_params%stencil_sizes

      global_domain_begin = block_params%domain_begin
      global_domain_end = block_params%domain_end
      global_domain_size = block_params%local_size

      ndims = comm_params%ndims

      start_dims = block_params%extended_block_begin_c+1
      dims = block_params%extended_block_dims

      call calculate_scaled_coefficients(ndims, block_params%dx, FDstencil_params)

      norm_array = [1e3,1e6,1e9,1e12]

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < solver_params%max_iter)
         norm_array(1) = 0.0
         do ii = begin(1), end(1)
            do jj = begin(2), end(2)
               index = [ii,jj]
               call IDX_XD(ndims, dims, index, global_index)

               old_val = block_params%matrix(global_index)

               call get_FD_coefficients_from_index(ndims, FDstencil_params%num_derivatives, FDstencil_params%stencil_sizes, &
                  start_dims, dims, index, FDstencil_params%scaled_stencil_coefficients, alphas, coefficients)

               dfxx => coefficients(1:FDstencil_params%num_stencil_elements)
               dfyy => coefficients(FDstencil_params%num_stencil_elements + 1:2 * FDstencil_params%num_stencil_elements)

               combined_stencils = dfxx + dfyy

               call f_analytical_poisson(ndims, funcs_params%rhs_func%output_size, start_dims, index, &
                  global_domain_begin, global_domain_end, global_domain_size, block_params%dx, f_val)

               call update_value_from_stencil(ndims, 1, 0, stencil_size, alphas, combined_stencils, &
                  dims, index, block_params%matrix, f_val(1), new_val)

               block_params%matrix(global_index) = new_val

               norm_array(1) = norm_array(1) + (abs(new_val - old_val))**2

            end do
         end do

         call sendrecv_data_neighbors(comm_params%comm, block_params, block_params%matrix_ptr)

         call check_convergence(comm_params%comm, solver_params%tol, solver_params%divergence_tol, &
            1.0/product(block_params%total_grid_size), norm_array, converged)
         if(converged == -1 .and. it > 0) then
            write(*,*) "Convergence failed"
            exit
         end if

         it = it + 1

      end do

      write(*,*) "Converged in ", it, " iterations"

   end subroutine GS_Method

   !!!! POISSON EQUATION FUNCTIONS !!!!

   pure subroutine initial_poisson(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      func_val = 20.0

   end subroutine initial_poisson

   pure subroutine boundary_poisson(ndims, num_elements, global_begin_indices, local_indices, &
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
         call u_analytical_poisson(ndims, num_elements, global_begin_indices, local_indices, &
            global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      end if

   end subroutine boundary_poisson

   pure subroutine f_analytical_poisson(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      real, dimension(ndims) :: point

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      func_val = -ndims*pi*pi*product(sin(pi*point))

      !func_val = 6.0*point(1) - 4 - 9*pi*pi*sin(3*pi*point(1))

   end subroutine f_analytical_poisson

   pure subroutine df_analytical_poisson(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      real, dimension(ndims) :: point

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      func_val = -ndims*pi*pi*pi*product(cos(pi*point))

   end subroutine df_analytical_poisson

   pure subroutine u_analytical_poisson(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      real, dimension(ndims) :: point

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      func_val = product(sin(pi*point))

      !func_val = point(1)*point(1)*point(1) - 2.0*point(1)*point(1) + sin(3*pi*point(1))

   end subroutine u_analytical_poisson

end module Poisson_module



