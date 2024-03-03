module solver_module
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM
   use utility_functions_module, only: IDX_XD, sleeper_function
   use block_module, only: block_type
   use rank_module, only: rank_type, communicate_step
   use finite_difference_module, only: FDstencil_type, apply_FDstencil
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper
   use functions_module, only : f_analytical_poisson_2d
   implicit none

   private
   public :: run_solver

contains

   function run_solver(parameters) result(result_array)
      type (rank_type), intent(inout) :: parameters
      integer :: iter, max_iter, converged

      real, dimension(4) :: result_array

      real :: local_norm, global_norm, previous_norm, relative_norm, max_tol, norm_scaling

      real, dimension(1) :: local_norm_array, global_norm_array

      norm_scaling = 1.0/product(parameters%grid_size) ! Scale depending on the number of grid points

      local_norm = 10e6
      global_norm = 10e9
      previous_norm = 10e18
      relative_norm = 10e36

      converged = 0
      max_iter = 1

      max_tol = 1.0e-6

      iter = 0
      do while (converged /= 1 .and. iter < max_iter)

         local_norm = Gauss_Seidel_iteration(parameters%ndims, parameters%block, parameters%FDstencil, &
            parameters%domain_begin)
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

         call communicate_step(parameters)

         iter = iter + 1
      end do

      result_array = [global_norm, relative_norm, REAL(converged,kind=8), REAL(iter,kind=8)]

   end function run_solver

   !> Gauss_Seidel_iteration with 2-norm
   !! Replace IDX_XD with a more efficient way to calculate the global index
   function Gauss_Seidel_iteration(ndims, block, FDstencil, global_domain_begin) result(norm)
      integer, intent(in) :: ndims
      type (block_type), intent(inout) :: block
      type (FDstencil_type), intent(in) :: FDstencil
      real, dimension(ndims) :: global_domain_begin

      integer, dimension(ndims) :: index, block_index
      integer :: ii, jj, local_index

      real :: stencil_val, f_val, new_val, norm

      stencil_val = 0.0
      f_val = 0.0
      new_val = 0.0
      norm = 0.0

      do ii = 2, block%size(1)-1
         do jj = 2, block%size(2)-1
            index = [ii,jj]
            block_index = block%begin + index - 1
            local_index = IDX_XD(ndims, block%size, index)

            stencil_val= apply_FDstencil(ndims, FDstencil, block, index)
            f_val = f_analytical_poisson_2d(ndims, global_domain_begin, block_index, FDstencil%dx) ! Could be an array instead of a function. Especially if we use jacobi-solvers. Depends on the size of the system.

            new_val = -stencil_val + f_val
            new_val = new_val / FDstencil%center_coefficient

            block%matrix(local_index) = new_val

            norm = norm + new_val**2
         end do
      end do

   end function Gauss_Seidel_iteration

end module solver_module
