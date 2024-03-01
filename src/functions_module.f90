module functions_module
use constants_module, only : pi
implicit none

private

abstract interface
   function myFuncInterface(ndims, global_domain_begin, global_indices, dx) result (func_val)
      integer, intent(in) :: ndims
      real, dimension(ndims), intent(in) :: global_domain_begin, dx
      integer, dimension(ndims), intent(in) :: global_indices

      real :: func_val

   end function myFuncInterface
end interface

public :: u_analytical_poisson_2d, f_analytical_poisson_2d

contains

function u_analytical_poisson_2d(ndims, global_domain_begin, global_indices, dx) result (func_val)
   integer, intent(in) :: ndims
   real, dimension(ndims), intent(in) :: global_domain_begin, dx
   integer, dimension(ndims), intent(in) :: global_indices

   real, dimension(ndims) :: point
   real :: func_val

   point = global_domain_begin + (global_indices(:) - 1) * dx(:)

   func_val = product(sin(pi*point))

end function u_analytical_poisson_2d

function f_analytical_poisson_2d(ndims, global_domain_begin, global_indices, dx) result (func_val)
   integer, intent(in) :: ndims
   real, dimension(ndims), intent(in) :: global_domain_begin, dx
   integer, dimension(ndims), intent(in) :: global_indices

   real, dimension(ndims) ::  point
   real :: func_val

   point = global_domain_begin + (global_indices(:) - 1) * dx(:)

   func_val = -2.0*pi*pi*product(sin(pi*point))

end function f_analytical_poisson_2d

end module functions_module
