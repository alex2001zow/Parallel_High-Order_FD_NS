module initilization_module
   use utility_functions_module
   implicit none

   private

   public :: initialize_block_2D

contains

   subroutine initialize_block_2D(ndim, global_dims, begin_block, dims, data, rank)
      integer, intent(in) :: ndim, rank
      integer, dimension(ndim), intent(in) :: global_dims, dims, begin_block
      real, dimension(product(dims)), intent(inout) :: data

      integer :: ii, jj, local_index, global_index

      do ii = 1, dims(1)
         do jj = 1, dims(2)
            local_index = IDX_XD(ndim,dims,[jj,ii])
            global_index = IDX_XD(ndim, global_dims, begin_block + [jj,ii] - 1)
            data(local_index) = test_function(global_index)
         end do
      end do

   end subroutine initialize_block_2D

   function test_function(input) result(output)
      integer, intent(in) :: input
      real :: output

      output = input

   end function test_function

end module initilization_module
