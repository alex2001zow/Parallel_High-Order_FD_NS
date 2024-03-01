module initialization_module
   use functions_module, only : u_analytical_poisson_2d
   use utility_functions_module, only : IDX_XD
   implicit none

   private

   public :: initialize_block_2D

contains

   !> Initialize the data array with the value of the function
   subroutine initialize_block_2D(ndims, global_dims, global_domain_begin, begin_block, dims, data)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: global_dims, dims, begin_block
      real, dimension(ndims), intent(in) :: global_domain_begin
      real, dimension(product(dims)), intent(inout) :: data

      integer :: ii, jj, local_index, global_index

      integer, dimension(ndims) :: index, block_index
      real, dimension(ndims) :: dx

      dx = 1.0 / (global_dims(:) - 1)

      do ii = 1, dims(1)
         do jj = 1, dims(2)
            index = [ii,jj]
            block_index = begin_block + index - 1
            local_index = IDX_XD(ndims, dims, index)
            global_index = IDX_XD(ndims, global_dims, block_index)
            data(local_index) = u_analytical_poisson_2d(ndims, global_domain_begin, block_index, dx)
         end do
      end do

      ! do ii = 2, dims(1)-1
      !    do jj = 2, dims(2)-1
      !       index = [ii,jj]
      !       local_index = IDX_XD(ndims, dims, index)
      !       data(local_index) = 0.0
      !    end do
      ! end do

   end subroutine initialize_block_2D


end module initialization_module
