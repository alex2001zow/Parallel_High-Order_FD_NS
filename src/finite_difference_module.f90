module finite_difference_module
   use utility_functions_module, only : IDX_XD
   implicit none

   private
   public :: finite_difference_stencil


contains

   !> Calculate the finite difference stencil coefficients.
   !! sum_{n=-alpha}^{beta} a_n * f(x + n*h)
   subroutine finite_difference_stencil(derivative_number, h, alpha, beta, stencil_coefficients)
      real, intent(in) :: h
      integer, intent(in) :: derivative_number, alpha, beta
      real, dimension((beta + alpha + 1)), intent(out) :: stencil_coefficients

      integer :: i, s, index, global_index
      real, dimension(alpha + beta + 1) :: taylor_coefficients
      real, dimension((alpha + beta + 1)**2) :: coefficient_matrix

      s = alpha + beta + 1

      ! Initialize the stencil coefficients. It starts as the rhs of the linear system.
      stencil_coefficients(:) = 0.0
      stencil_coefficients(derivative_number) = 1.0

      index = 1
      do i = alpha, beta
         global_index = IDX_XD(2, [s,s], [1, index])
         call calculate_taylor_coefficients(i, h, s, taylor_coefficients)
         coefficient_matrix(global_index:global_index + s - 1) = taylor_coefficients
         index = index + 1
      end do

      ! Solve the linear system to get the stencil coefficients. Use LAPACK. Save the output to stencil_coefficients. We need to make sure the stride is correct!
      ! coefficient_matrix * x = stencil_coefficients

   end subroutine finite_difference_stencil

   !> Calculate the Taylor coefficients up to order s
   !! sum_{i=0}^{s-1} (n*h)^i / i!
   subroutine calculate_taylor_coefficients(n, h, s, coefficients)
      integer, intent(in) :: n, s
      real, intent(in) :: h
      real, dimension(s), intent(out) :: coefficients

      integer :: i
      real :: numerator, denominator

      numerator = 1.0
      denominator = 1.0
      do i = 1, s
         coefficients(i) = numerator / denominator
         numerator = numerator * n*h
         denominator = denominator * i
      end do

   end subroutine calculate_taylor_coefficients

end module finite_difference_module
