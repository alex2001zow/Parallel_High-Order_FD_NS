module solver_module
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM

   use utility_functions_module, only: IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers
   use initialization_module, only: write_function_to_block

   use comm_module, only: comm_type
   use block_module, only: block_type, sendrecv_data_neighbors
   use FD_module, only: FDstencil_type, apply_FDstencil, calculate_scaled_coefficients, get_FD_coefficients_from_index
   use functions_module, only: FunctionPair, FunctionPtrType

   use mpi_wrapper_module, only: all_reduce_mpi_wrapper
   implicit none

   private

   abstract interface
      pure subroutine SystemSolverInterface(ndims, num_elements, stencil_size, alphas, num_derivatives, &
         stencil_coefficients, combined_stencils, dims, index, matrix, f_val, residual, Jac_residual)
         integer, intent(in) :: ndims, num_elements, num_derivatives
         integer, dimension(ndims), intent(in) :: stencil_size, alphas, dims, index
         real, dimension(:), target, intent(inout) :: stencil_coefficients
         real, dimension(:), intent(in) :: matrix
         real, dimension(:), intent(in) :: f_val

         real, dimension(product(stencil_size)), intent(inout) :: combined_stencils
         real, dimension(num_elements), intent(inout) :: residual, Jac_residual ! The approximate residual function and Jacobian of the residual for a given point.

      end subroutine SystemSolverInterface
   end interface

   type SolverPtrType
      procedure(SystemSolverInterface), pointer, nopass :: solver => null()
   end type SolverPtrType

   type SolverParamsType
      real :: tol, divergence_tol
      integer :: max_iter, solver_type
   end type SolverParamsType

   type ResultType
      integer :: converged
      integer :: iterations

      real :: global_u_norm
      real :: global_r_norm

      real :: relative_u_norm
      real :: relative_r_norm
   end type ResultType

   enum, bind(C)
      enumerator :: JSolver = 1
      enumerator :: GSSolver = 2
   end enum

   public :: SolverPtrType, set_SystemSolver_pointer
   public :: SolverParamsType, set_SolverParamsType
   public :: ResultType, set_ResultType, print_ResultType
   public :: choose_iterative_solver
   public :: check_convergence

contains

   !> Set the pointer to the system solver
   pure subroutine set_SystemSolver_pointer(system_func, SystemSolver)
      procedure(SystemSolverInterface) :: system_func
      type(SolverPtrType), intent(inout) :: SystemSolver

      SystemSolver%solver => system_func

   end subroutine set_SystemSolver_pointer

   !> Set the solver parameters. Solver type is Jacobi=1 and Gauss-Seidel=2
   pure subroutine set_SolverParamsType(tol, divergence_tol, max_iter, solver_type, SolverParams)
      real, intent(in) :: tol, divergence_tol
      integer, intent(in) :: max_iter, solver_type
      type(SolverParamsType), intent(out) :: SolverParams

      SolverParams%tol = tol
      SolverParams%divergence_tol = divergence_tol
      SolverParams%max_iter = max_iter
      SolverParams%solver_type = solver_type

   end subroutine set_SolverParamsType

   !> Set the result type
   pure subroutine set_ResultType(converged, iterations, global_u_norm, global_r_norm, relative_u_norm, relative_r_norm, result)
      integer, intent(in) :: converged, iterations
      real, intent(in) :: global_u_norm, global_r_norm, relative_u_norm, relative_r_norm
      type(ResultType), intent(out) :: result

      result%converged = converged
      result%iterations = iterations
      result%global_u_norm = global_u_norm
      result%global_r_norm = global_r_norm
      result%relative_u_norm = relative_u_norm
      result%relative_r_norm = relative_r_norm

   end subroutine set_ResultType

   !> Print the result type
   subroutine print_ResultType(Result)
      type(ResultType), intent(in) :: Result

      write(*,"(A, I12.1)") "Converged: ", Result%converged
      write(*,"(A, I12.1)") "Iterations: ", Result%iterations
      write(*,"(A, E12.6)") "Global u_norm: ", Result%global_u_norm
      write(*,"(A, E12.6)") "Global r_norm: ", Result%global_r_norm
      write(*,"(A, E12.6)") "Relative u_norm: ", Result%relative_u_norm
      write(*,"(A, E12.6)") "Relative r_norm: ", Result%relative_r_norm

   end subroutine print_ResultType

   subroutine choose_iterative_solver(comm_in, block_in, FDstencil_in, &
      functions_in, SystemSolver_in, SolverParams_in, Result_out)
      type(comm_type), target, intent(inout) :: comm_in
      type(block_type), target, intent(inout) :: block_in
      type(FDstencil_type), target, intent(inout) :: FDstencil_in
      type(FunctionPair), target, intent(in) :: functions_in
      type(SolverPtrType), target, intent(in) :: SystemSolver_in
      type(SolverParamsType), intent(in) :: SolverParams_in

      type(ResultType), intent(inout) :: Result_out

      call calculate_scaled_coefficients(comm_in%ndims, block_in%extended_grid_dx, FDstencil_in)

      select case (SolverParams_in%solver_type)
       case(GSSolver)
         call run_GSS_solver(comm_in, block_in, FDstencil_in, functions_in, SystemSolver_in, SolverParams_in, Result_out)
       case(JSolver)
         print *, "Jacobi solver not implemented!"
         stop
         !call run_JS_solver(comm_in, block_in, FDstencil_in, functions_in, SystemSolver_in, SolverParams_in, &
         !   begin, end, Result_out)
       case default
         print *, "Solver not implemented or does not exist!"
         stop
      end select

   end subroutine choose_iterative_solver

   !> Run the Gauss-Seidel solver.
   subroutine run_GSS_solver(comm_in, block_in, FDstencil_in, functions_in, SystemSolver_in, SolverParams_in, Result_out)
      type(comm_type), target, intent(inout) :: comm_in
      type(block_type), target, intent(inout) :: block_in
      type(FDstencil_type), target, intent(inout) :: FDstencil_in
      type(FunctionPair), target, intent(in) :: functions_in
      type(SolverPtrType), target, intent(in) :: SystemSolver_in
      type(SolverParamsType), intent(in) :: SolverParams_in

      type(ResultType), intent(inout) :: Result_out

      integer :: iter, converged
      real :: norm_scaling
      real, dimension(10) :: norm_array

      ! Scale depending on the number of grid points
      norm_scaling = 1.0/product(block_in%extended_grid_size)

      ! local_norm, global_norm, previous_norm, relative_norm
      norm_array = 1e3

      converged = 0
      iter = 0
      do while (converged /= 1 .and. iter < SolverParams_in%max_iter)

         call GS_iteration(block_in%ndims, block_in%num_elements, &
            FDstencil_in, block_in%extended_block_begin_c+1, block_in%extended_block_dims, &
            block_in%global_begin_c+1, block_in%matrix, functions_in, block_in%extended_grid_dx, SystemSolver_in, &
            block_in%domain_begin, block_in%domain_end, block_in%extended_grid_size, &
            block_in%block_begin_c+1, block_in%block_end_c, norm_array(1))

         call check_convergence(comm_in%comm, SolverParams_in%tol, SolverParams_in%divergence_tol, iter, &
            norm_array, converged)
         if(converged == -1) then
            exit
         end if

         call sendrecv_data_neighbors(comm_in%comm, block_in, block_in%matrix_ptr)

         iter = iter + 1
      end do

      !call set_ResultType(converged, iter, norm_array(2), norm_array(4), Result_out)

   end subroutine run_GSS_solver

   !> Gauss Seidel iteration with 2-norm
   subroutine GS_iteration(ndims, num_elements, FDstencil, &
      start_dims, dims, global_begin, matrix, functions, dx, SystemSolver, &
      global_domain_begin, global_domain_end, global_domain_size, begin, end, norm)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: start_dims, dims, global_begin, global_domain_size, begin, end
      real, dimension(ndims), intent(in) :: dx

      type(FDstencil_type), intent(inout) :: FDstencil
      type(FunctionPair), intent(in) :: functions
      type(SolverPtrType), intent(in) :: SystemSolver
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end

      real, dimension(:), intent(inout) :: matrix
      real, intent(inout) :: norm

      integer, dimension(ndims) :: local_dims, index, alphas, betas
      integer :: global_index, local_index, start_index, end_index
      real, dimension(functions%rhs_func%output_size) :: f_val, old_val, new_val, residual, Jac_residual

      real, dimension(FDstencil%num_stencil_elements) :: combined_stencil_coefficients

      real, contiguous, dimension(:), pointer :: pointer_to_coefficients

      local_dims = end - begin + 1

      norm = 0.0

      do global_index = 1, product(local_dims)

         call IDX_XD_INV(ndims, local_dims, global_index, index)
         index = begin + index - 1

         call IDX_XD(ndims, dims, index, local_index)
         start_index = (local_index - 1) * functions%rhs_func%output_size + 1
         end_index = local_index * functions%rhs_func%output_size
         old_val = matrix(start_index:end_index)

         ! Calculate the function
         call functions%rhs_func%func(ndims, functions%rhs_func%output_size, global_begin, index, &
            global_domain_begin, global_domain_end, global_domain_size, dx, f_val)

         call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
            start_dims, dims, index, FDstencil%scaled_stencil_coefficients, alphas, betas, pointer_to_coefficients)

         ! Call the system solver
         call SystemSolver%solver(ndims, functions%rhs_func%output_size, &
            FDstencil%stencil_sizes, alphas, FDstencil%num_derivatives, &
            pointer_to_coefficients, &
            combined_stencil_coefficients, dims, index, matrix, f_val, residual, Jac_residual)

         call Newtons_iteration(old_val, residual, Jac_residual, new_val)

         matrix(start_index:end_index) = new_val

         norm = norm + sum(abs(residual)**2) ! 2-norm
      end do

   end subroutine GS_iteration

   !> Newtons method iteration
   elemental subroutine Newtons_iteration(x_old, residual, Jac_residual, x_new)
      real, intent(in) :: x_old, residual, Jac_residual
      real, intent(out) :: x_new

      x_new = x_old - residual / (Jac_residual + 1e-6)

   end subroutine Newtons_iteration

   ! Calculate the relative difference
   elemental subroutine calculate_relative_difference(x_new, x_old, relative_difference)
      real, intent(in) :: x_new, x_old
      real, intent(out) :: relative_difference

      relative_difference = abs(x_new - x_old) / (abs(x_new) + 1e-6)

   end subroutine calculate_relative_difference

   !> Check for convergence
   !! norm_array(1) = local u_diff_norm
   !! norm_array(2) = local r_diff_norm
   !! norm_array(3) = local u_norm
   !! norm_array(4) = local r_norm
   !! norm_array(5) = global u_diff_norm
   !! norm_array(6) = global r_diff_norm
   !! norm_array(7) = global u_norm
   !! norm_array(8) = global r_norm
   !! norm_array(9) = previous relative u_norm
   !! norm_array(10) = previous relative r_norm
   subroutine check_convergence(comm, tol, divergence_tol, it, norm_array, converged)
      integer, intent(in) :: comm, it
      real, intent(in) ::  tol, divergence_tol
      real, dimension(10), intent(inout) :: norm_array
      integer, intent(out) :: converged

      real :: relative_u_norm, relative_r_norm

      ! Find the global norm from all blocks
      call all_reduce_mpi_wrapper(norm_array(1:4), norm_array(5:8), 4, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm)

      ! Calculate ||u_new - u_old||/||u_new||
      relative_u_norm = norm_array(5)/(norm_array(7)+1e-6)

      ! Calculate ||r_new - r_old||/||r_new||
      relative_r_norm = norm_array(6)/(norm_array(8)+1e-6)

      if(relative_u_norm < tol .and. relative_r_norm < tol) then
         converged = 1
      else if(relative_u_norm > norm_array(9) + divergence_tol .or. &
         relative_r_norm > norm_array(10) + divergence_tol .and. it > 5) then
         converged = -1
      else
         converged = 0 ! Still converging
      end if

      norm_array(9) = relative_u_norm ! Overwrite the previous relative u_norm
      norm_array(10) = relative_r_norm ! Overwrite the previous relative r_norm
   end subroutine

end module solver_module
