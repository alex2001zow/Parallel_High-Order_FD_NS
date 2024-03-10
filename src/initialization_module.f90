module initialization_module
   use functions_module, only : u_analytical_poisson_2d
   use utility_functions_module, only : IDX_XD
   implicit none

   private

   public :: initialize_block_2D

contains

   !> Initialize the buffer array with the value of the function
   subroutine initialize_block_2D(ndims, global_dims, global_domain_begin, begin_block, dims, buffer, dx, rank)
      integer, intent(in) :: ndims, rank
      integer, dimension(ndims), intent(in) :: global_dims, dims, begin_block
      real, dimension(ndims), intent(in) :: global_domain_begin, dx
      real, dimension(product(dims)), intent(inout) :: buffer

      integer :: ii, jj, local_index

      integer, dimension(ndims) :: index, block_index

      do ii = 1, dims(1)
         do jj = 1, dims(2)
            index = [ii,jj]
            block_index = begin_block + index - 1
            local_index = IDX_XD(ndims, dims, index)
            !buffer(local_index) = rank!IDX_XD(ndims, global_dims, block_index)
            buffer(local_index) = u_analytical_poisson_2d(ndims, global_domain_begin, block_index, dx)
         end do
      end do

      ! do ii = 2, dims(1)-1
      !    do jj = 2, dims(2)-1
      !       index = [ii,jj]
      !       local_index = IDX_XD(ndims, dims, index)
      !       buffer(local_index) = 1.0
      !    end do
      ! end do

   end subroutine initialize_block_2D


end module initialization_module
