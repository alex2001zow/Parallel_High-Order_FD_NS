module utility_functions_module
   implicit none

   private
   public :: IDX_XD, print_cartesian_grid

contains

   !> Global index from ndims, dims, indices.
   !! 1D: ndims = 1, dims = [N], indices = [i], global_index = i
   !! 2D: ndims = 2, dims = [N, M], indices = [i, j], global_index = i + N * (j - 1)
   !! 3D: ndims = 3, dims = [N, M, K], indices = [i, j, k], global_index = i + N * (j - 1) + N * M * (k - 1)
   !! 4D: ndims = 4, dims = [N, M, K, L], indices = [i, j, k, l], global_index = i + N * (j - 1) + N * M * (k - 1) + N * M * K * (l - 1)
   !! etc.
   function IDX_XD(ndims, dims, indices) result(global_index)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims, indices
      integer :: global_index
      integer :: d, product

      global_index = indices(1) - 1
      product = 1

      do d = 2, ndims
         product = product * dims(d - 1)
         global_index = global_index + (indices(d) - 1) * product
      end do

      ! Adjust index for 1-based indexing of Fortran
      global_index = global_index + 1
   end function IDX_XD

   !> A routine to print the cartesian grid. Just for debugging and understanding
   subroutine print_cartesian_grid(ndim, pn)
      implicit none
      integer, intent(in) :: ndim  ! Number of dimensions
      integer, dimension(ndim), intent(in) :: pn  ! Processors in each dimension
      integer :: i, j, k, idx
      integer, dimension(ndim) :: indices  ! Array to hold the indices for IDX_XD

      print *, "Cartesian processor grid with dimension:", ndim

      select case(ndim)
       case(1)
         ! 1D Grid
         do i = 1, pn(1)
            indices(1) = i
            idx = IDX_XD(ndim, pn, indices) - 1
            write(*, '(I4, 1X)', advance="no") idx
         end do
         print *

       case(2)
         ! 2D Grid
         do j = 1, pn(2)
            do i = 1, pn(1)
               indices = [i, j]  ! Set current indices
               idx = IDX_XD(ndim, pn, indices) - 1
               write(*, '(I4, 1X)', advance="no") idx
            end do
            print *  ! New line for the next row
         end do

       case(3)
         ! 3D Grid (printed as slices of 2D grids)
         do k = 1, pn(3)
            print *, "Slice", k, ":"
            do j = 1, pn(2)
               do i = 1, pn(1)
                  indices = [i, j, k]  ! Set current indices
                  idx = IDX_XD(ndim, pn, indices) - 1
                  write(*, '(I4, 1X)', advance="no") idx
               end do
               print *  ! New line for the next row
            end do
            if (k < pn(3)) then
               print *  ! Separate slices with a blank line for readability
            endif
         end do

      end select

   end subroutine print_cartesian_grid

end module utility_functions_module
