module solver_module
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM

   use utility_functions_module, only: IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers

   use comm_module, only: comm_type
   use block_module, only: block_type, sendrecv_data_neighbors
   use FD_module, only: FDstencil_type, apply_FDstencil, calculate_scaled_coefficients, get_FD_coefficients_from_index

   use mpi_wrapper_module, only: all_reduce_mpi_wrapper
   implicit none

   private

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

   public :: SolverParamsType, set_SolverParamsType
   public :: ResultType, set_ResultType, print_ResultType
   public :: check_convergence

contains

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
