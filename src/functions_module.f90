module functions_module
implicit none

private

abstract interface
   pure subroutine FunctionInterface(ndims, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, intent(inout) :: func_val

   end subroutine FunctionInterface
end interface

type :: FunctionPtrType
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

pure subroutine set_function_pointers(initial_func, boundary_func, rhs_func, rhs_derivative_func, analytical_func, funcPair)
   procedure(FunctionInterface) :: initial_func, boundary_func, rhs_func, rhs_derivative_func, analytical_func
   type(FunctionPair), intent(inout) :: funcPair

   funcPair%initial_condition_func%func => initial_func
   funcPair%boundary_condition_func%func => boundary_func
   funcPair%rhs_func%func => rhs_func
   funcPair%rhs_derivative_func%func => rhs_derivative_func
   funcPair%analytical_func%func => analytical_func

end subroutine set_function_pointers

pure subroutine calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)
   integer, intent(in) :: ndims
   integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices
   real, dimension(ndims), intent(in) :: global_domain_begin, dx
   real, dimension(ndims), intent(inout) :: point

   integer, dimension(ndims) :: block_index

   ! Calculate where in global space the point is
   block_index = global_begin_indices + local_indices - 1

   ! Calculate the point in the global domain
   point = global_domain_begin + (block_index - 1) * dx

end subroutine calculate_point

end module functions_module
