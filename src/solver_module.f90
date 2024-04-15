module solver_module
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM

   use utility_functions_module, only: IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers
   use initialization_module, only: write_function_to_block

   use rank_module, only: rank_type, communicate_step
   use comm_module, only: comm_type
   use block_module, only: block_type
   use FD_module, only: FDstencil_type, apply_FDstencil, determine_alpha, alpha_2_global, global_2_start_end
   use functions_module, only: FunctionPair, FunctionPtrType

   use mpi_wrapper_module, only: all_reduce_mpi_wrapper
   implicit none

   private

   abstract interface
      pure subroutine SystemSolverInterface(ndims, stencil_size, alphas, num_derivatives, stencil_coefficients, &
         dims, index, matrix, f_val, df_val, F, J)
         integer, intent(in) :: ndims, num_derivatives
         integer, dimension(ndims), intent(in) :: stencil_size, alphas, dims, index
         real, dimension(product(stencil_size)), intent(in) :: stencil_coefficients
         real, dimension(product(dims)), intent(in) :: matrix
         real, intent(in) :: f_val, df_val

         real, intent(inout) :: F, J ! The approximate function and Jacobian value for a given point.

      end subroutine SystemSolverInterface
   end interface

   type SolverPtrType
      procedure(SystemSolverInterface), pointer, nopass :: solver => null()
   end type SolverPtrType

   enum, bind(C)
      enumerator :: DefaultSolver = 0
      enumerator :: GSSolver = 1
      enumerator :: JSolver = 2
      enumerator :: OtherSolver = 3
   end enum

   public :: SolverPtrType, set_SystemSolver_pointer, run_solver

contains

   pure subroutine set_SystemSolver_pointer(system_func, SystemSolver)
      procedure(SystemSolverInterface) :: system_func
      type(SolverPtrType), intent(inout) :: SystemSolver

      SystemSolver%solver => system_func

   end subroutine set_SystemSolver_pointer

   !> Run the solver. Split this into a function for each solver type?
   subroutine run_solver(parameters_in, comm_in, block_in, FDstencil_in, functions_in, SystemSolver_in, begin, end, &
      result_array)
      type(rank_type), target, intent(in) :: parameters_in
      type(comm_type), target, intent(inout) :: comm_in
      type(block_type), target, intent(inout) :: block_in
      type(FDstencil_type), target, intent(in) :: FDstencil_in
      type(FunctionPair), target, intent(in) :: functions_in
      type(SolverPtrType), target, intent(in) :: SystemSolver_in
      integer, dimension(parameters_in%ndims), intent(in) :: begin, end

      real, dimension(4), intent(inout) :: result_array

      integer :: iter, max_iter, converged
      real, dimension(:), pointer :: ptr_temp_array, ptr_matrix
      real :: local_norm, global_norm, previous_norm, relative_norm, max_tol, norm_scaling
      real, dimension(1) :: local_norm_array, global_norm_array
      integer :: enum_solver

      enum_solver = GSSolver

      select case(enum_solver)
       case(GSSolver)
         ! Do nothing
       case(JSolver)
         ! Copy the matrix to the temp array
         block_in%temp_array = block_in%matrix
         ! Copy the matrix to the f_array (just to create a matrix of the same size as the block to read from)
         block_in%f_array = block_in%temp_array
         ! Write the function to the f_array
         call write_function_to_block(parameters_in%ndims, parameters_in%domain_begin, parameters_in%domain_end, &
            parameters_in%grid_size, block_in%begin, block_in%size, block_in%f_array, &
            FDstencil_in%dx, functions_in%rhs_func)
       case default
         ! Do nothing
      end select

      ! Initialize pointers
      ptr_matrix => block_in%matrix
      ptr_temp_array => block_in%temp_array

      norm_scaling = 1.0/product(parameters_in%grid_size) ! Scale depending on the number of grid points

      local_norm = 1e3
      global_norm = 1e6
      previous_norm = 1e9
      relative_norm = 1e18

      converged = 0
      max_iter = 10000

      max_tol = 1e-6

      iter = 0
      do while (converged /= 1 .and. iter < max_iter)

         ! Perhaps we should make a seperate function depending on the solver. So we do not have to check the solver type in the loop
         select case(enum_solver)
          case(GSSolver)
            call GS_iteration(parameters_in%ndims, FDstencil_in, block_in%size*0 + 1, block_in%size, &
               block_in%begin, ptr_matrix, functions_in, SystemSolver_in, &
               parameters_in%domain_begin, parameters_in%domain_end, parameters_in%grid_size, begin, end, local_norm)
          case(JSolver)
            ! We set the matrix as the temp array to make sure we have the newest data in the matrix and not the temp_array
            call Jacobi_iteration(parameters_in%ndims, FDstencil_in, block_in%size, &
               block_in%begin, ptr_temp_array, block_in%f_array, ptr_matrix, local_norm)
          case default
            print *, "Solver not implemented or does not exist!"
            stop
         end select

         local_norm_array(1) = local_norm

         call all_reduce_mpi_wrapper(local_norm_array, global_norm_array, 1, &
            int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_in%comm)

         global_norm = global_norm_array(1) * norm_scaling
         global_norm = sqrt(global_norm)
         ! Check for divergence: if the norm increases
         if (global_norm > previous_norm .and. iter > 0) then
            converged = -1  ! Indicate divergence
            exit  ! Exit the loop
         end if

         relative_norm = abs((global_norm - previous_norm) / (previous_norm + 1e-6))
         if (previous_norm > 0.0 .and. relative_norm < max_tol) then
            converged = 1
         end if

         previous_norm = global_norm

         call communicate_step(parameters_in%ndims, comm_in, block_in, ptr_matrix)

         select case(enum_solver)
          case(GSSolver)
            ! Do nothing
          case(JSolver)
            ! Swap the pointers using a function
            call swap_pointers(ptr_matrix, ptr_temp_array)
          case default
            print *, "Solver not implemented or does not exist!"
            stop
         end select

         iter = iter + 1
      end do

      result_array = [global_norm, relative_norm, real(converged,kind=8), real(iter,kind=8)]

   end subroutine run_solver

   !> Jacobi iteration with 2-norm. Should be like gauss-seidel but with a copy of the matrix. Fix it when Gauss-Seidel works
   subroutine Jacobi_iteration(ndims, FDstencil, dims, global_begin, matrix, f_array, temp_array, norm)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims, global_begin
      real, dimension(product(dims)), intent(inout) :: matrix, f_array, temp_array
      type(FDstencil_type), intent(in) :: FDstencil

      real, intent(inout) :: norm

      integer, dimension(ndims) :: begin, end, local_dims, index, block_index, alphas
      integer :: global_index, local_index

      real :: stencil_val, f_val, old_val, new_val

      ! This is just for this stencil. This should be dynamic depending on the stencil
      begin = 2
      end = dims - 1
      local_dims = end - begin + 1

      norm = 0.0

      !$omp parallel do reduction(+:norm) private(stencil_val, f_val, old_val, new_val, &
      !$omp& index, local_index, block_index, alphas) shared(ndims, FDstencil, dims, global_begin, &
      !$omp& matrix, f_array, temp_array, begin, end, local_dims) default(none)
      do global_index = 1, product(local_dims)
         call IDX_XD_INV(ndims, local_dims, global_index, index)
         index = begin + index - 1
         call IDX_XD(ndims, dims, index, local_index)
         block_index = global_begin + index - 1

         call determine_alpha(ndims, FDstencil%stencil_sizes, begin, end, index, alphas)

         call apply_FDstencil(ndims, FDstencil%stencil_sizes, alphas, FDstencil%stencil_coefficients, &
            dims, index, matrix, stencil_val)
         f_val = f_array(local_index)

         new_val = stencil_val - f_val
         old_val = matrix(local_index)

         temp_array(local_index) = new_val

         norm = norm + (new_val-old_val)*(new_val-old_val)
      end do

   end subroutine Jacobi_iteration

   !> Gauss Seidel iteration with 2-norm
   pure subroutine GS_iteration(ndims, FDstencil, start_dims, dims, global_begin, matrix, functions, SystemSolver, &
      global_domain_begin, global_domain_end, global_domain_size, begin, end, norm)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: start_dims, dims, global_begin, global_domain_size, begin, end

      type(FDstencil_type), intent(in) :: FDstencil
      type(FunctionPair), intent(in) :: functions
      type(SolverPtrType), intent(in) :: SystemSolver
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end

      real, dimension(product(dims)), intent(inout) :: matrix
      real, intent(inout) :: norm

      integer, dimension(ndims) :: local_dims, index, block_index, alphas
      real, dimension(ndims) :: point
      integer :: global_index, local_index, coefficient_global_index, coefficients_start_index, coefficients_end_index
      real :: f_val, df_val, old_val, F, J, new_val

      local_dims = end - begin + 1

      norm = 0.0

      do global_index = 1, product(local_dims)

         call IDX_XD_INV(ndims, local_dims, global_index, index)
         index = begin + index - 1
         block_index = global_begin + index - 1

         call IDX_XD(ndims, dims, index, local_index)
         old_val = matrix(local_index)

         ! Calculate the point in the global domain
         point = global_domain_begin + (block_index - 1) * FDstencil%dx

         ! Determine the alphas for the current index
         call determine_alpha(ndims, FDstencil%stencil_sizes, start_dims, dims, index, alphas)

         ! Find the coefficients from alpha.
         call alpha_2_global(ndims, FDstencil%stencil_sizes, alphas, coefficient_global_index)
         call global_2_start_end(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, coefficient_global_index, &
            coefficients_start_index, coefficients_end_index)

         ! Calculate the function and derivative values
         call functions%rhs_func%func(ndims, global_domain_begin, global_domain_end, global_domain_size, &
            block_index, FDstencil%dx, point, f_val)
         call functions%rhs_derivative_func%func(ndims, global_domain_begin, global_domain_end, global_domain_size, &
            block_index, FDstencil%dx, point, df_val)

         ! Call the system solver
         call SystemSolver%solver(ndims, FDstencil%stencil_sizes, alphas, FDstencil%num_derivatives, &
            FDstencil%stencil_coefficients(coefficients_start_index:coefficients_end_index), &
            dims, index, matrix, f_val, df_val, F, J)

         call Newtons_iteration(old_val, F, J, new_val)

         matrix(local_index) = new_val

         norm = norm + (new_val-old_val)**2
      end do

   end subroutine GS_iteration

   !> Newtons method iteration
   elemental subroutine Newtons_iteration(x, F, J, x_new)
      real, intent(in) :: x, F, J
      real, intent(inout) :: x_new

      x_new = x - F / (J + 1e-6) ! J + 1e-6 to avoid division by zero.

   end subroutine Newtons_iteration

end module solver_module
