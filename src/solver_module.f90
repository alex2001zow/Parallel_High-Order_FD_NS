module solver_module
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM
   use functions_module, only: FunctionPtrType
   use utility_functions_module, only: IDX_XD, sleeper_function, swap_pointers
   use initialization_module, only: write_function_to_block
   use block_module, only: block_type
   use rank_module, only: rank_type, communicate_step
   use finite_difference_module, only: FDstencil_type, apply_FDstencil
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

   function run_solver(parameters) result(result_array)
      type (rank_type), target, intent(inout) :: parameters
      integer :: iter, max_iter, converged

      real, dimension(4) :: result_array
      real, dimension(:), pointer :: ptr_temp_array, ptr_matrix

      real :: local_norm, global_norm, previous_norm, relative_norm, max_tol, norm_scaling

      real, dimension(1) :: local_norm_array, global_norm_array

      integer :: enum_solver

      enum_solver = JSolver

      select case(enum_solver)
       case(GSSolver)
         ! Do nothing
       case(JSolver)
         ! Copy the matrix to the temp array
         parameters%block%temp_array = parameters%block%matrix
         ! Copy the matrix to the f_array (just to create a matrix of the same size as the block to read from)
         parameters%block%f_array = parameters%block%temp_array
         ! Write the function to the f_array
         call write_function_to_block(parameters%ndims, parameters%domain_begin, &
            parameters%block%begin, parameters%block%size, parameters%block%f_array, &
            parameters%FDstencil%dx, parameters%funcs%rhs_func)
       case default
         ! Do nothing
      end select

      ! Initialize pointers
      ptr_matrix => parameters%block%matrix
      ptr_temp_array => parameters%block%temp_array

      norm_scaling = 1.0/product(parameters%grid_size - 1) ! Scale depending on the number of grid points

      local_norm = 10e6
      global_norm = 10e9
      previous_norm = 10e18
      relative_norm = 10e36

      converged = 0
      max_iter = 10000

      max_tol = 1.0e-6

      iter = 0
      do while (converged /= 1 .and. iter < max_iter)

         ! Perhaps we should make a seperate function depending on the solver. So we do not have to check the solver type in the loop
         select case(enum_solver)
          case(GSSolver)
            local_norm = Gauss_Seidel_iteration(parameters%ndims, parameters%FDstencil,parameters%block%size, &
               parameters%block%begin, parameters%block%matrix, parameters%funcs%rhs_func, parameters%domain_begin)
          case(JSolver)
            ! We set the matrix as the temp array to make sure we have the newest data in the matrix and not the temp_array
            local_norm = Jacobi_iteration(parameters%ndims, parameters%FDstencil,parameters%block%size, &
               parameters%block%begin, ptr_temp_array, parameters%block%f_array, ptr_matrix)
          case default
            print *, "Solver not implemented or does not exist!"
            stop
         end select

         local_norm_array(1) = local_norm

         call all_reduce_mpi_wrapper(local_norm_array, global_norm_array, 1, &
            INT(MPI_DOUBLE_PRECISION,kind=8), INT(MPI_SUM,kind=8), parameters%comm%comm)

         global_norm = global_norm_array(1) * norm_scaling
         global_norm = sqrt(global_norm) ! Not sure if needed

         relative_norm = abs((global_norm - previous_norm) / previous_norm)

         if (previous_norm > 0.0 .and. relative_norm < max_tol) then
            converged = 1
         end if

         previous_norm = global_norm

         call communicate_step(parameters%ndims, parameters%comm, parameters%block, ptr_matrix)

         ! Swap the pointers using a function
         call swap_pointers(ptr_matrix, ptr_temp_array)

         iter = iter + 1
      end do

      result_array = [global_norm, relative_norm, REAL(converged,kind=8), REAL(iter,kind=8)]

   end function run_solver

   !! Replace IDX_XD with a more efficient way to calculate the global index depending on the number of dimensions. This is true for both Gauss_Seidel_iteration and Jacobi_iteration

   !> Jacobi iteration with 2-norm
   function Jacobi_iteration(ndims, FDstencil, dims, begin, matrix, f_array, temp_array) result(norm)
      integer, intent(in) :: ndims
      integer, dimension(ndims) :: dims, begin
      real, dimension(product(dims)), intent(inout) :: matrix, f_array, temp_array
      type (FDstencil_type), intent(in) :: FDstencil

      integer, dimension(ndims) :: index, block_index
      integer :: ii, jj, local_index

      real :: stencil_val, f_val, new_val, norm

      stencil_val = 0.0
      f_val = 0.0
      new_val = 0.0
      norm = 0.0

      !$omp parallel do reduction(+:norm) private(stencil_val, f_val, new_val, &
      !$omp& index, block_index, local_index) shared(ndims, FDstencil, dims, begin, &
      !$omp& matrix, f_array, temp_array) default(none) collapse(2)
      do ii = 2, dims(1)-1
         do jj = 2, dims(2)-1
            index = [ii,jj]
            block_index = begin + index - 1
            local_index = IDX_XD(ndims, dims, index)

            stencil_val= apply_FDstencil(ndims, FDstencil, dims, matrix, index)
            f_val = f_array(local_index)

            new_val = -stencil_val + f_val
            new_val = new_val / FDstencil%center_coefficient

            temp_array(local_index) = new_val

            norm = norm + new_val**2
         end do
      end do

   end function Jacobi_iteration

   !> Gauss_Seidel_iteration with 2-norm
   function Gauss_Seidel_iteration(ndims, FDstencil, dims, begin, matrix, f_func, global_domain_begin) result(norm)
      integer, intent(in) :: ndims
      integer, dimension(ndims) :: dims, begin
      real, dimension(product(dims)), intent(inout) :: matrix
      type (FDstencil_type), intent(in) :: FDstencil
      type(FunctionPtrType), intent(in) :: f_func
      real, dimension(ndims) :: global_domain_begin

      integer, dimension(ndims) :: index, block_index
      integer :: ii, jj, local_index

      real :: stencil_val, f_val, new_val, norm

      stencil_val = 0.0
      f_val = 0.0
      new_val = 0.0
      norm = 0.0

      do ii = 2, dims(1)-1
         do jj = 2, dims(2)-1
            index = [ii,jj]
            block_index = begin + index - 1
            local_index = IDX_XD(ndims, dims, index)

            stencil_val= apply_FDstencil(ndims, FDstencil, dims, matrix, index)
            f_val = f_func%func(ndims, global_domain_begin, block_index, FDstencil%dx)

            new_val = -stencil_val + f_val
            new_val = new_val / FDstencil%center_coefficient

            matrix(local_index) = new_val

            norm = norm + new_val**2
         end do
      end do

   end function Gauss_Seidel_iteration

end module solver_module
