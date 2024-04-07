module functions_module
implicit none

private

abstract interface
   pure subroutine FunctionInterface(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
      integer, intent(in) :: ndims
      real, dimension(ndims), intent(in) :: domain_start, domain_end, dx, point
      integer, dimension(ndims), intent(in) :: domain_size, domain_indices

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
   type(FunctionPtrType) :: analytical_func
end type FunctionPair

public :: FunctionPtrType, FunctionPair, set_function_pointers

contains

pure subroutine set_function_pointers(initial_func, boundary_func, rhs_func, analytical_func, funcPair)
   !type(FunctionPtrType), intent(in) :: initial_func, boundary_func, rhs_func, analytical_func
   procedure(FunctionInterface) :: initial_func, boundary_func, rhs_func, analytical_func
   type(FunctionPair), intent(out) :: funcPair

   ! Initialize the function pointers to null
   funcPair%initial_condition_func%func => initial_func
   funcPair%boundary_condition_func%func => boundary_func
   funcPair%rhs_func%func => rhs_func
   funcPair%analytical_func%func => analytical_func

end subroutine set_function_pointers

end module functions_module
