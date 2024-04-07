module solver_module
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM

   use utility_functions_module, only: IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers
   use initialization_module, only: write_function_to_block

   use rank_module, only: rank_type, communicate_step
   use comm_module, only: comm_type
   use block_module, only: block_type
   use FD_module, only: FDstencil_type, apply_FDstencil, determine_alpha
   use functions_module, only: FunctionPair, FunctionPtrType

   use mpi_wrapper_module, only: all_reduce_mpi_wrapper
   implicit none

   private

   enum, bind(C)
      enumerator :: DefaultSolver = 0
      enumerator :: GSSolver = 1
      enumerator :: JSolver = 2
      enumerator :: OtherSolver = 3
   end enum

   public :: run_solver

contains

   !> Run the solver. Split this into a function for each solver type?
   subroutine run_solver(parameters_in, comm_in, block_in, FDstencil_in, functions_in, result_array)
      type(rank_type), target, intent(in) :: parameters_in
      type(comm_type), target, intent(inout) :: comm_in
      type(block_type), target, intent(inout) :: block_in
      type(FDstencil_type), target, intent(in) :: FDstencil_in
      type(FunctionPair), target, intent(in) :: functions_in

      integer :: iter, max_iter, converged

      real, dimension(4), intent(out) :: result_array
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

      local_norm = 10e3
      global_norm = 10e6
      previous_norm = 10e9
      relative_norm = 10e18

      converged = 0
      max_iter = 100000

      max_tol = 1.0e-6

      iter = 0
      do while (converged /= 1 .and. iter < max_iter)

         ! Perhaps we should make a seperate function depending on the solver. So we do not have to check the solver type in the loop
         select case(enum_solver)
          case(GSSolver)
            local_norm = Gauss_Seidel_iteration(parameters_in%ndims, FDstencil_in, block_in%size, &
               block_in%begin, ptr_matrix, functions_in%rhs_func, parameters_in%domain_begin, &
               parameters_in%domain_end, parameters_in%grid_size)
          case(JSolver)
            ! We set the matrix as the temp array to make sure we have the newest data in the matrix and not the temp_array
            local_norm = Jacobi_iteration(parameters_in%ndims, FDstencil_in, block_in%size, &
               block_in%begin, ptr_temp_array, block_in%f_array, ptr_matrix)
          case default
            print *, "Solver not implemented or does not exist!"
            stop
         end select

         local_norm_array(1) = local_norm

         call all_reduce_mpi_wrapper(local_norm_array, global_norm_array, 1, &
            int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_in%comm)

         global_norm = global_norm_array(1) * norm_scaling
         global_norm = sqrt(global_norm) ! Not sure if needed

         relative_norm = abs((global_norm - previous_norm) / previous_norm)

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

   !> Jacobi iteration with 2-norm. Should be like gauss-seidel but with a copy of the matrix. Fix it when gauss-Seidel works
   function Jacobi_iteration(ndims, FDstencil, dims, global_begin, matrix, f_array, temp_array) result(norm)
      integer, intent(in) :: ndims
      integer, dimension(ndims) :: dims, global_begin
      real, dimension(product(dims)), intent(inout) :: matrix, f_array, temp_array
      type(FDstencil_type), intent(in) :: FDstencil

      integer, dimension(ndims) :: begin, end, local_dims, index, block_index, alphas
      integer :: global_index, local_index

      real :: stencil_val, f_val, old_val, new_val, norm

      ! This is just for this stencil. This should be dynamic depending on the stencil
      begin = 2
      end = dims - 1
      local_dims = end - begin + 1

      norm = 0.0

      !$omp parallel do reduction(+:norm) private(stencil_val, f_val, old_val, new_val, &
      !$omp& index, local_index, block_index, alphas) shared(ndims, FDstencil, dims, global_begin, &
      !$omp& matrix, f_array, temp_array, begin, end, local_dims) default(none)
      do global_index = 1, product(local_dims)
         index = IDX_XD_INV(ndims, local_dims, global_index)
         index = begin + index - 1
         local_index = IDX_XD(ndims, dims, index)
         block_index = global_begin + index - 1

         alphas = determine_alpha(ndims, FDstencil%stencil_sizes, begin, end, index)

         stencil_val = apply_FDstencil(ndims, FDstencil%stencil_sizes, alphas, FDstencil%stencil_coefficients, &
            dims, index, matrix)
         f_val = f_array(local_index)

         new_val = stencil_val - f_val
         old_val = matrix(local_index)

         temp_array(local_index) = new_val

         norm = norm + (new_val-old_val)*(new_val-old_val)
      end do

   end function Jacobi_iteration

   !> Gauss_Seidel_iteration with 2-norm
   function Gauss_Seidel_iteration(ndims, FDstencil, dims, global_begin, matrix, f_func, &
      global_domain_begin, global_domain_end, global_domain_size) result(norm)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims, global_begin, global_domain_size

      type(FDstencil_type), intent(in) :: FDstencil
      type(FunctionPtrType), intent(in) :: f_func
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end

      real, dimension(product(dims)), intent(inout) :: matrix

      integer, dimension(ndims) :: begin, end, local_dims, index, block_index, alphas
      real, dimension(ndims) :: point
      integer :: global_index, local_index, coefficient_index, coefficients_start_index, coefficients_end_index
      real :: stencil_val, f_val, old_val, new_val, norm, omega

      begin = 2
      end = dims - 1
      local_dims = end - begin + 1

      norm = 0.0

      omega = 1e-2

      do global_index = 1, product(local_dims)
         index = IDX_XD_INV(ndims, local_dims, global_index)
         index = begin + index - 1
         local_index = IDX_XD(ndims, dims, index)
         block_index = global_begin + index - 1

         point = global_domain_begin + (block_index - 1) * FDstencil%dx

         alphas = determine_alpha(ndims, FDstencil%stencil_sizes, begin, dims, index)

         coefficient_index = IDX_XD(ndims, FDstencil%stencil_sizes, FDstencil%stencil_sizes-alphas-1)
         coefficients_start_index = coefficient_index*FDstencil%num_derivatives*FDstencil%num_stencil_elements + 1
         coefficients_end_index = (coefficient_index + 1)*FDstencil%num_derivatives*FDstencil%num_stencil_elements

         stencil_val = apply_FDstencil(ndims, FDstencil%stencil_sizes, alphas, &
            FDstencil%stencil_coefficients(coefficients_start_index:coefficients_end_index), dims, index, matrix)
         call f_func%func(ndims, global_domain_begin, global_domain_end, global_domain_size, &
            block_index, FDstencil%dx, point, f_val)

         old_val = matrix(local_index)
         new_val = stencil_val - f_val

         matrix(local_index) = (1.0-omega)*old_val + omega*new_val

         norm = norm + (new_val-old_val)*(new_val-old_val)
      end do

   end function Gauss_Seidel_iteration

end module solver_module
