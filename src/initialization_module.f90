module initialization_module
   use functions_module, only: FunctionPtrType, FunctionPair
   use utility_functions_module, only: IDX_XD
   implicit none

   private

   public :: write_function_to_block, write_initial_condition_and_boundary

contains

   !> Write function values to a block
   subroutine write_function_to_block(ndims, global_domain_begin, begin_block, dims, buffer, dx, func)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims, begin_block
      real, dimension(ndims), intent(in) :: global_domain_begin, dx
      real, dimension(product(dims)), intent(inout) :: buffer
      type(FunctionPtrType), intent(in) :: func

      integer :: ii, jj, local_index
      integer, dimension(ndims) :: index, block_index

      do ii = 1, dims(1)
         do jj = 1, dims(2)
            index = [ii,jj]
            block_index = begin_block + index - 1
            local_index = IDX_XD(ndims, dims, index)
            ! Write function value to buffer
            buffer(local_index) = func%func(ndims, global_domain_begin, block_index, dx)
         end do
      end do

   end subroutine write_function_to_block

   !> Initialize the buffer with the inital condition and the boundary condition
   subroutine write_initial_condition_and_boundary(ndims, global_domain_begin, begin_block, dims, buffer, dx, &
      inital_func, boundary_func)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims, begin_block
      real, dimension(ndims), intent(in) :: global_domain_begin, dx
      real, dimension(product(dims)), intent(inout) :: buffer
      type(FunctionPtrType), intent(in) :: inital_func, boundary_func

      integer :: ii, jj, local_index
      integer, dimension(ndims) :: index, block_index

      do ii = 1, dims(1)
         do jj = 1, dims(2)
            index = [ii,jj]
            block_index = begin_block + index - 1
            local_index = IDX_XD(ndims, dims, index)

            ! Write initial condition
            buffer(local_index) = inital_func%func(ndims, global_domain_begin, block_index, dx)
            ! Overwrite with boundary conditions

            buffer(local_index) = boundary_func%func(ndims, global_domain_begin, block_index, dx)
         end do
      end do

   end subroutine write_initial_condition_and_boundary

end module initialization_module
