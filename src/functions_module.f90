module functions_module
implicit none

private

public :: calculate_point

contains

pure subroutine calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)
   integer, intent(in) :: ndims
   integer, dimension(:), intent(in) :: global_begin_indices, local_indices
   real, dimension(:), intent(in) :: global_domain_begin, dx
   real, dimension(:), intent(out) :: point

   integer, dimension(ndims) :: block_index

   ! Calculate where in global space the point is
   block_index = global_begin_indices + local_indices

   ! Calculate the point in the global domain
   point = global_domain_begin + (block_index - 1) * dx

end subroutine calculate_point

end module functions_module
