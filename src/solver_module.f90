module solver_module
   use rank_parameters_module, only: rank_type, communicate_step
   implicit none


   private
   public :: run_solver

contains

   function run_solver(parameters) result(converged)
      type (rank_type), intent(inout) :: parameters
      integer :: iter, max_iter
      logical :: converged

      converged = .false.
      max_iter = 1000

      iter = 0
      do while (converged .neqv. .true. .and. iter < max_iter)
         call communicate_step(parameters)
         converged = .true.
         iter = iter + 1
      end do

   end function


end module solver_module
