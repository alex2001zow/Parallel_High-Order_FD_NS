module functions_module
use constants_module, only : pi
implicit none

private

enum, bind(C)
   enumerator :: DefaultModel = 0
   enumerator :: PoissonModel = 1
   enumerator :: AnotherModel = 2
   ! Add more models as needed
end enum

abstract interface
   subroutine FunctionInterface(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
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

public :: PoissonModel
public :: FunctionPtrType, FunctionPair, set_function_pointers

contains

subroutine set_function_pointers(enum_val, funcPair)
   integer, intent(in) :: enum_val
   type(FunctionPair), intent(out) :: funcPair

   ! Initialize the function pointers to null
   funcPair%initial_condition_func%func => null()
   funcPair%boundary_condition_func%func => null()
   funcPair%rhs_func%func => null()
   funcPair%analytical_func%func => null()

   ! Check for the matching name and set the function pointers
   select case(enum_val)
    case(PoissonModel)
      funcPair%initial_condition_func%func => initial_poisson
      funcPair%boundary_condition_func%func => boundary_poisson
      funcPair%rhs_func%func => f_analytical_poisson
      funcPair%analytical_func%func => u_analytical_poisson
    case default
      write(*,*) 'No matching function pair for the given name.'
   end select
end subroutine set_function_pointers

subroutine initial_poisson(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
   integer, intent(in) :: ndims
   real, dimension(ndims), intent(in) :: domain_start, domain_end, dx, point
   integer, dimension(ndims), intent(in) :: domain_size, domain_indices

   real, intent(inout) :: func_val

   func_val = 20.0

end subroutine initial_poisson

subroutine boundary_poisson(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
   integer, intent(in) :: ndims
   real, dimension(ndims), intent(in) :: domain_start, domain_end, dx, point
   integer, dimension(ndims), intent(in) :: domain_size, domain_indices

   real, intent(inout) :: func_val

   integer :: i
   logical :: at_boundary
   at_boundary = .false.

   do i = 1, ndims
      if (domain_indices(i) == 1 .or. domain_indices(i) == domain_size(i)) then
         at_boundary = .true.
         exit
      end if
   end do

   if (at_boundary) then
      call u_analytical_poisson(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
   end if

end subroutine boundary_poisson

subroutine u_analytical_poisson(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
   integer, intent(in) :: ndims
   real, dimension(ndims), intent(in) :: domain_start, domain_end, dx, point
   integer, dimension(ndims), intent(in) :: domain_size, domain_indices

   real, intent(inout) :: func_val

   func_val = product(sin(pi*point))

end subroutine u_analytical_poisson

subroutine f_analytical_poisson(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
   integer, intent(in) :: ndims
   real, dimension(ndims), intent(in) :: domain_start, domain_end, dx, point
   integer, dimension(ndims), intent(in) :: domain_size, domain_indices

   real, intent(inout) :: func_val

   func_val = -ndims*pi*pi*product(sin(pi*point))

end subroutine f_analytical_poisson

end module functions_module
