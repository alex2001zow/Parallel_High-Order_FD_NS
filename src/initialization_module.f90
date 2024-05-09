module initialization_module
   use functions_module, only: FunctionPtrType, FunctionPair, calculate_point
   use utility_functions_module, only: IDX_XD_INV
   implicit none

   private

   public :: write_function_to_block, write_initial_condition_and_boundary

contains

   !> Write function values to a block
   pure subroutine write_function_to_block(ndims, num_elements, global_domain_begin, global_domain_end, global_domain_size, &
      global_begin_block, dims, buffer, dx, func)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) ::  global_domain_size, global_begin_block, dims
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx
      real, dimension(product(dims)*num_elements), intent(inout) :: buffer
      type(FunctionPtrType), intent(in) :: func

      integer :: global_index, index_start, index_end
      integer, dimension(ndims) :: index
      real, dimension(num_elements) :: func_val

      do global_index = 1, product(dims)
         call IDX_XD_INV(ndims, dims, global_index, index)
         index_start = (global_index - 1) * num_elements + 1
         index_end = global_index * num_elements

         ! Write function value to buffer
         call func%func(ndims, num_elements, global_begin_block, index, &
            global_domain_begin, global_domain_end, global_domain_size, dx, func_val)

         buffer(index_start:index_end) = func_val

      end do

   end subroutine write_function_to_block

   !> Initialize the buffer with the inital condition and the boundary condition
   pure subroutine write_initial_condition_and_boundary(ndims, num_elements, &
      global_domain_begin, global_domain_end, global_domain_size, &
      global_begin_block, dims, buffer, dx, functions_in)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) ::  global_domain_size, global_begin_block, dims
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx
      real, dimension(product(dims)*num_elements), intent(inout) :: buffer
      type(FunctionPair), intent(in) :: functions_in

      integer :: global_index, index_start, index_end
      integer, dimension(ndims) :: index
      real, dimension(num_elements) :: func_val

      do global_index = 1, product(dims)
         call IDX_XD_INV(ndims, dims, global_index, index)
         index_start = (global_index - 1) * num_elements + 1
         index_end = global_index * num_elements

         ! Write initial condition
         call functions_in%initial_condition_func%func(ndims, num_elements, global_begin_block, index, &
            global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
         buffer(index_start:index_end) = func_val

         ! Overwrite with boundary conditions
         call functions_in%boundary_condition_func%func(ndims, num_elements, global_begin_block, index, &
            global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
         buffer(index_start:index_end) = func_val
      end do

   end subroutine write_initial_condition_and_boundary

end module initialization_module
