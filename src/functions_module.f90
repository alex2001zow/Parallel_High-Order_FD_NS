module functions_module
use constants_module, only : pi
implicit none

private

enum, bind(C)
   enumerator :: DefaultModel = 0
   enumerator :: Poisson2D = 1
   enumerator :: AnotherModel = 2
   ! Add more models as needed
end enum

abstract interface
   function FunctionInterface(ndims, global_domain_begin, global_indices, dx) result (func_val)
      integer, intent(in) :: ndims
      real, intent(in) :: global_domain_begin(ndims), dx(ndims)
      integer, intent(in) :: global_indices(ndims)

      real :: func_val

   end function FunctionInterface
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

public :: Poisson2D
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
    case(Poisson2D)
      funcPair%initial_condition_func%func => u_analytical_poisson_2d
      funcPair%boundary_condition_func%func => u_analytical_poisson_2d
      funcPair%rhs_func%func => f_analytical_poisson_2d
      funcPair%analytical_func%func => u_analytical_poisson_2d
    case default
      write(*,*) 'No matching function pair for the given name.'
   end select
end subroutine set_function_pointers

function initial_poisson_2d(ndims, global_domain_begin, global_indices, dx) result (func_val)
   integer, intent(in) :: ndims
   real, dimension(ndims), intent(in) :: global_domain_begin, dx
   integer, dimension(ndims), intent(in) :: global_indices

   real, dimension(ndims) :: point
   real :: func_val

   point = global_domain_begin + (global_indices(:) - 1) * dx(:)

   func_val = 1.0

end function initial_poisson_2d

function boundary_poisson_2d(ndims, global_domain_begin, global_indices, dx) result (func_val)
   integer, intent(in) :: ndims
   real, dimension(ndims), intent(in) :: global_domain_begin, dx
   integer, dimension(ndims), intent(in) :: global_indices

   real, dimension(ndims) :: point
   real :: func_val

   point = global_domain_begin + (global_indices(:) - 1) * dx(:)

   func_val = product(sin(pi*point))

end function boundary_poisson_2d

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
