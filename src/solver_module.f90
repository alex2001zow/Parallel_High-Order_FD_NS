module solver_module
   use rank_parameters_module
   use initilization_module
   implicit none


   private
   public :: run_solver

contains

   function run_solver(parameters) result(converged)
      type (rank_struct), intent(inout) :: parameters
      logical :: converged

      converged = .false.

   end function


end module solver_module
