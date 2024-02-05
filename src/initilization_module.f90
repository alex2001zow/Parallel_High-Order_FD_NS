module initilization_module
   use utility_functions_module
   implicit none

   private

   public :: test_block_structure_2D

contains

   subroutine test_block_structure_2D(ndim, dims, data, rank)
      integer, intent(in) :: ndim, rank
      integer, dimension(ndim), intent(in) :: dims
      real, dimension(product(dims)), intent(inout) :: data

      integer :: ii, jj, global_index

      do ii = 1, dims(1)
         do jj = 1, dims(2)
            global_index = IDX_XD(ndim,dims,[jj,ii])
            data(global_index) = rank
         end do
      end do

   end subroutine test_block_structure_2D

end module initilization_module
