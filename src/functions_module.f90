module functions_module
use constants_module, only : pi
implicit none

private
public :: u_analytical_poisson_2d, f_analytical_poisson_2d


contains

function u_analytical_poisson_2d(ndims, grid_size, indices) result (u)
   integer, intent(in) :: ndims
   integer, dimension(ndims), intent(in) :: grid_size, indices

   real, dimension(ndims) :: dx, point
   real :: u

   dx(:) = 1.0/REAL(grid_size(:),kind=8)

   point = (indices(:) - 1) * dx(:)

   u = product(sin(pi*point))

end function u_analytical_poisson_2d

function f_analytical_poisson_2d(ndims, grid_size, indices) result (f)
   integer, intent(in) :: ndims
   integer, dimension(ndims), intent(in) :: grid_size, indices

   real, dimension(ndims) :: dx, point
   real :: f

   dx(:) = 1.0/REAL(grid_size(:),kind=8)

   point = (indices(:) - 1) * dx(:)

   f = -2.0*pi*pi*product(sin(pi*point))

end function f_analytical_poisson_2d

end module functions_module
