module initialization_module
   use functions_module, only : u_analytical_poisson_2d
   use utility_functions_module, only : IDX_XD
   implicit none

   private

   public :: initialize_block_2D

contains

   !> Initialize the data array with the value of the function
   subroutine initialize_block_2D(ndims, global_dims, begin_block, dims, data, rank)
      integer, intent(in) :: ndims, rank
      integer, dimension(ndims), intent(in) :: global_dims, dims, begin_block
      real, dimension(product(dims)), intent(inout) :: data

      integer :: ii, jj, local_index, global_index

      do ii = 1, dims(1)
         do jj = 1, dims(2)
            local_index = IDX_XD(ndims, dims, [ii,jj])
            global_index = IDX_XD(ndims, global_dims, begin_block + [ii,jj] - 1)
            data(local_index) = global_index!u_analytical_poisson_2d(ndims, global_dims, begin_block + [ii,jj])
         end do
      end do

   end subroutine initialize_block_2D


end module initialization_module
