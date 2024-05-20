!> This module contains implementations of the necessary functions for the multigrid method
module multigrid_module
   use utility_functions_module, only: IDX_XD, IDX_XD_INV
   implicit none

contains

   !> Simple injection restriction operator for 2D problems
   subroutine injection_restriction_2D(ndims, fine_dims, coarse_dims, fine_matrix, coarse_matrix)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: fine_dims, coarse_dims
      real, dimension(:), intent(in) :: fine_matrix
      real, dimension(:), intent(inout) :: coarse_matrix

      integer :: ii, jj, ii_fine, jj_fine, coarse_index, fine_index
      integer, dimension(ndims) :: coarse_indices, fine_indices


      ! Make it parallel
      do ii = 1, coarse_dims(1)
         ii_fine = 2*ii - 1
         do jj = 1, coarse_dims(2)
            jj_fine = 2*jj - 1
            coarse_indices = [ii,jj]
            fine_indices = [ii_fine,jj_fine]
            call IDX_XD(ndims,coarse_dims,coarse_indices,coarse_index)
            call IDX_XD(ndims,fine_dims,fine_indices,fine_index)
            coarse_matrix(coarse_index) = fine_matrix(fine_index)
         end do
      end do

   end subroutine injection_restriction_2D

   !> Simple linear interpolation operator for 2D problems
   subroutine linear_prolongation_2D(ndims, fine_dims, coarse_dims, fine_matrix, coarse_matrix)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: fine_dims, coarse_dims
      real, dimension(:), intent(in) :: coarse_matrix
      real, dimension(:), intent(inout) :: fine_matrix

      integer :: ii, jj, ii_fine, jj_fine, coarse_index, fine_index
      integer, dimension(ndims) :: coarse_indices, fine_indices
      real :: weight

      ! Make it parallel
      do ii = 1, coarse_dims(1)
         ii_fine = 2*ii - 1
         do jj = 1, coarse_dims(2)
            jj_fine = 2*jj - 1
            coarse_indices = [ii,jj]
            fine_indices = [ii_fine,jj_fine]
            call IDX_XD(ndims,coarse_dims,coarse_indices,coarse_index)
            call IDX_XD(ndims,fine_dims,fine_indices,fine_index)
            fine_matrix(fine_index) = coarse_matrix(coarse_index)
            ! Interpolate in the x direction
            if (ii < coarse_dims(1)) then
               fine_indices = [ii_fine+1,jj_fine]
               call IDX_XD(ndims,fine_dims,fine_indices,fine_index)
               weight = 0.5
               fine_matrix(fine_index) = fine_matrix(fine_index) + weight*coarse_matrix(coarse_index)
            end if
            ! Interpolate in the y direction
            if (jj < coarse_dims(2)) then
               fine_indices = [ii_fine,jj_fine+1]
               call IDX_XD(ndims,fine_dims,fine_indices,fine_index)
               weight = 0.5
               fine_matrix(fine_index) = fine_matrix(fine_index) + weight*coarse_matrix(coarse_index)
            end if
            ! Interpolate in the x and y directions
            if (ii < coarse_dims(1) .and. jj < coarse_dims(2)) then
               fine_indices = [ii_fine+1,jj_fine+1]
               call IDX_XD(ndims,fine_dims,fine_indices,fine_index)
               weight = 0.25
               fine_matrix(fine_index) = fine_matrix(fine_index) + weight*coarse_matrix(coarse_index)
            end if
         end do
      end do

   end subroutine linear_prolongation_2D

end module multigrid_module
