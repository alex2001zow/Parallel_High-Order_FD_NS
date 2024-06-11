module functions_module
implicit none

private

abstract interface
   subroutine FunctionInterface(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

   end subroutine FunctionInterface
end interface

type :: FunctionPtrType
   integer :: output_size
   procedure(FunctionInterface), pointer, nopass :: func => null()
end type FunctionPtrType

type :: FunctionPair
   type(FunctionPtrType) :: initial_condition_func
   type(FunctionPtrType) :: boundary_condition_func
   type(FunctionPtrType) :: rhs_func
   type(FunctionPtrType) :: rhs_derivative_func
   type(FunctionPtrType) :: analytical_func
end type FunctionPair

public :: FunctionPtrType, FunctionPair, set_function_pointers
public :: calculate_point

contains

subroutine set_function_pointers(output_size, initial_func, boundary_func, &
   rhs_func, rhs_derivative_func, analytical_func, funcPair)
   integer, dimension(5), intent(in) :: output_size
   procedure(FunctionInterface) :: initial_func, boundary_func, rhs_func, rhs_derivative_func, analytical_func
   type(FunctionPair), intent(out) :: funcPair

   funcPair%initial_condition_func%output_size = output_size(1)
   funcPair%initial_condition_func%func => initial_func

   funcPair%boundary_condition_func%output_size = output_size(2)
   funcPair%boundary_condition_func%func => boundary_func

   funcPair%rhs_func%output_size = output_size(3)
   funcPair%rhs_func%func => rhs_func

   funcPair%rhs_derivative_func%output_size = output_size(4)
   funcPair%rhs_derivative_func%func => rhs_derivative_func

   funcPair%analytical_func%output_size = output_size(5)
   funcPair%analytical_func%func => analytical_func

end subroutine set_function_pointers

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
