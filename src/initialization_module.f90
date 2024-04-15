module initialization_module
   use functions_module, only: FunctionPtrType, FunctionPair
   use utility_functions_module, only: IDX_XD_INV
   implicit none

   private

   public :: write_function_to_block, write_initial_condition_and_boundary

contains

   !> Write function values to a block
   pure subroutine write_function_to_block(ndims, global_domain_begin, global_domain_end, global_domain_size, &
      begin_block, dims, buffer, dx, func)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) ::  global_domain_size, begin_block, dims
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx
      real, dimension(product(dims)), intent(inout) :: buffer
      type(FunctionPtrType), intent(in) :: func

      integer :: global_index
      integer, dimension(ndims) :: index, block_index
      real, dimension(ndims) :: point
      real :: func_val

      do global_index = 1, product(dims)
         call IDX_XD_INV(ndims, dims, global_index, index)
         block_index = begin_block + index - 1

         point = global_domain_begin + (block_index - 1) * dx

         ! Write function value to buffer
         call func%func(ndims, global_domain_begin, global_domain_end, global_domain_size, block_index, dx, point, func_val)
         buffer(global_index) = func_val

      end do

   end subroutine write_function_to_block

   !> Initialize the buffer with the inital condition and the boundary condition
   pure subroutine write_initial_condition_and_boundary(ndims, global_domain_begin, global_domain_end, global_domain_size, &
      begin_block, dims, buffer, dx, inital_func, boundary_func)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) ::  global_domain_size, begin_block, dims
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx
      real, dimension(product(dims)), intent(inout) :: buffer
      type(FunctionPtrType), intent(in) :: inital_func, boundary_func

      integer :: global_index
      integer, dimension(ndims) :: index, block_index
      real, dimension(ndims) :: point
      real :: func_val

      do global_index = 1, product(dims)
         call IDX_XD_INV(ndims, dims, global_index, index)
         block_index = begin_block + index - 1

         point = global_domain_begin + (block_index - 1) * dx

         ! Write initial condition
         call inital_func%func(ndims, global_domain_begin, global_domain_end, global_domain_size, &
            block_index, dx, point, func_val)
         buffer(global_index) = func_val

         ! Overwrite with boundary conditions
         call boundary_func%func(ndims, global_domain_begin, global_domain_end, global_domain_size, &
            block_index, dx, point, func_val)
         buffer(global_index) = func_val
      end do

   end subroutine write_initial_condition_and_boundary

end module initialization_module
