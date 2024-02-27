module solver_module
   use mpi, only: MPI_REAL, MPI_SUM
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

      integer, dimension(2) :: result_array

      real :: local_norm, global_norm, max_tol

      real, dimension(1) :: local_norm_array, global_norm_array

      local_norm = 999.0
      global_norm = 0.0

      converged = 0
      max_iter = 1

      max_tol = 1.0e-6

      iter = 0
      do while (converged /= 1 .and. iter < max_iter)

         local_norm = Gauss_Seidel_iteration(parameters%ndims, parameters%block, parameters%FDstencil, parameters%grid_size)

         local_norm_array(1) = local_norm

         !call all_reduce_mpi_wrapper(local_norm_array, global_norm_array, 1, &
         !   INT(MPI_REAL,kind=8), INT(MPI_SUM,kind=8), parameters%comm%comm)

         global_norm = global_norm_array(1)

         !global_norm = sqrt(global_norm / product(parameters%grid_size))

         if(global_norm < max_tol) then
            !converged = 1
         end if

         !call communicate_step(parameters)

         iter = iter + 1
      end do

      result_array = [converged, iter]

   end function run_solver

   !> Gauss_Seidel_iteration with 2-norm
   !! Replace IDX_XD with a more efficient way to calculate the global index
   function Gauss_Seidel_iteration(ndims, block, FDstencil, grid_size) result(norm)
      integer, intent(in) :: ndims
      type (block_type), intent(inout) :: block
      type (FDstencil_type), intent(in) :: FDstencil
      integer, dimension(ndims), intent(in) :: grid_size

      integer, dimension(ndims) ::index
      integer :: ii, jj, global_index

      real :: f_val, norm, val

      f_val = 0.0
      norm = 0.0
      val = 0.0

      do ii = 2, block%size(1)-1
         do jj = 2, block%size(2)-1
            index = [ii,jj]
            global_index = IDX_XD(ndims, block%size, index)

            val = apply_FDstencil(ndims, FDstencil, block, index)/product(grid_size**ndims)
            f_val = f_analytical_poisson_2d(ndims, grid_size, block%begin + [ii,jj])

            val = val + f_val
            norm = norm + val**2
            block%matrix(global_index) = val
         end do
      end do

   end function Gauss_Seidel_iteration

end module solver_module
