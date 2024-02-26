module solver_module
   use utility_functions_module, only: IDX_XD_REVERSED
   use block_module, only: block_type
   use rank_module, only: rank_type, communicate_step
   use finite_difference_module, only: FDstencil_type, apply_FDstencil
   implicit none

   private
   public :: run_solver

contains

   function run_solver(parameters) result(converged)
      type (rank_type), intent(inout) :: parameters
      integer :: iter, max_iter
      logical :: converged

      real :: norm, norm_new, norm_reduce

      norm = 0.0
      norm_new = 0.0
      norm_reduce = 0.0

      converged = .false.
      max_iter = 101

      iter = 0
      do while (converged .neqv. .true. .and. iter < max_iter)

         !norm_new = Gauss_Seidel_iteration(parameters%ndims, parameters%block, parameters%FDstencil)

         !! MPI allreduce norm to check for convergence

         call communicate_step(parameters)

         !converged = .true.
         iter = iter + 1
      end do

   end function

   !> Gauss_Seidel_iteration with 2-norm
   !! Replace IDX_XD with a more efficient way to calculate the global index
   function Gauss_Seidel_iteration(ndims, block, FDstencil) result(norm)
      integer, intent(in) :: ndims
      type (block_type), intent(inout) :: block
      type (FDstencil_type), intent(in) :: FDstencil

      integer :: ii, jj, global_index

      real :: f, norm, val

      f = 0.0
      norm = 0.0
      val = 0.0

      do ii = 2, block%size(1)-1
         do jj = 2, block%size(2)-1
            global_index = IDX_XD_REVERSED(ndims, block%size, [ii,jj])
            val = apply_FDstencil(ndims, FDstencil, block, [ii,jj])
            norm = norm + val**2
            block%matrix(global_index) = val
         end do
      end do

   end function Gauss_Seidel_iteration

end module solver_module
