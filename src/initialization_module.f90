module initialization_module
   use functions_module, only: FunctionPtrType, FunctionPair, calculate_point
   use utility_functions_module, only: IDX_XD_INV
   implicit none

   private

   public :: write_function_to_block, write_initial_condition_and_boundary

contains

   !> Write function values to a block
   pure subroutine write_function_to_block(ndims, global_domain_begin, global_domain_end, global_domain_size, &
      global_begin_block, dims, buffer, dx, func)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) ::  global_domain_size, global_begin_block, dims
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx
      real, dimension(product(dims)), intent(inout) :: buffer
      type(FunctionPtrType), intent(in) :: func

      integer :: global_index
      integer, dimension(ndims) :: index
      real :: func_val

      do global_index = 1, product(dims)
         call IDX_XD_INV(ndims, dims, global_index, index)

         ! Write function value to buffer
         call func%func(ndims, global_begin_block, index, &
            global_domain_begin, global_domain_end, global_domain_size, dx, func_val)

         buffer(global_index) = func_val

      end do

   end subroutine write_function_to_block

   !> Initialize the buffer with the inital condition and the boundary condition
   pure subroutine write_initial_condition_and_boundary(ndims, global_domain_begin, global_domain_end, global_domain_size, &
      global_begin_block, dims, buffer, dx, functions_in)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) ::  global_domain_size, global_begin_block, dims
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx
      real, dimension(product(dims)), intent(inout) :: buffer
      type(FunctionPair), intent(in) :: functions_in

      integer :: global_index
      integer, dimension(ndims) :: index
      real :: func_val

      do global_index = 1, product(dims)
         call IDX_XD_INV(ndims, dims, global_index, index)

         ! Write initial condition
         call functions_in%initial_condition_func%func(ndims, global_begin_block, index, &
            global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
         buffer(global_index) = func_val

         ! Overwrite with boundary conditions
         call functions_in%boundary_condition_func%func(ndims, global_begin_block, index, &
            global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
         buffer(global_index) = func_val
      end do

   end subroutine write_initial_condition_and_boundary

end module initialization_module
