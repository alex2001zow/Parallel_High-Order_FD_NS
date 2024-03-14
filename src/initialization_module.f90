module initialization_module
   use functions_module, only: FunctionPtrType, FunctionPair
   use utility_functions_module, only: IDX_XD
   implicit none

   private

   public :: write_function_to_block, write_initial_condition_and_boundary

contains

   !> Write function values to a block
   subroutine write_function_to_block(ndims, global_domain_begin, global_domain_end, global_domain_size, &
      begin_block, dims, buffer, dx, func)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) ::  global_domain_size, begin_block, dims
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx
      real, dimension(product(dims)), intent(inout) :: buffer
      type(FunctionPtrType), intent(in) :: func

      integer :: ii, jj, local_index
      integer, dimension(ndims) :: index, block_index
      real, dimension(ndims) :: point
      real :: func_val

      do ii = 1, dims(1)
         do jj = 1, dims(2)
            index = [ii,jj]
            block_index = begin_block + index - 1
            local_index = IDX_XD(ndims, dims, index)

            point = global_domain_begin + (block_index - 1) * dx

            ! Write function value to buffer
            call func%func(ndims, global_domain_begin, global_domain_end, global_domain_size, block_index, dx, point, func_val)
            buffer(local_index) = func_val
         end do
      end do

   end subroutine write_function_to_block

   !> Initialize the buffer with the inital condition and the boundary condition
   subroutine write_initial_condition_and_boundary(ndims, global_domain_begin, global_domain_end, global_domain_size, &
      begin_block, dims, buffer, dx, inital_func, boundary_func)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) ::  global_domain_size, begin_block, dims
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx
      real, dimension(product(dims)), intent(inout) :: buffer
      type(FunctionPtrType), intent(in) :: inital_func, boundary_func

      integer :: ii, jj, local_index
      integer, dimension(ndims) :: index, block_index
      real, dimension(ndims) :: point
      real :: func_val

      do ii = 1, dims(1)
         do jj = 1, dims(2)
            index = [ii,jj]
            block_index = begin_block + index - 1
            local_index = IDX_XD(ndims, dims, index)

            point = global_domain_begin + (block_index - 1) * dx

            ! Write initial condition
            call inital_func%func(ndims, global_domain_begin, global_domain_end, global_domain_size, &
               block_index, dx, point, func_val)
            buffer(local_index) = func_val

            ! Overwrite with boundary conditions
            call boundary_func%func(ndims, global_domain_begin, global_domain_end, global_domain_size, &
               block_index, dx, point, func_val)
            buffer(local_index) = func_val
         end do
      end do

   end subroutine write_initial_condition_and_boundary

end module initialization_module
