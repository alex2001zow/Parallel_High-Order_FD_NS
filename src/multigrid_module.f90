!> This module contains implementations of the necessary functions for the multigrid method
module multigrid_module
   use utility_functions_module, only: IDX_XD, IDX_XD_INV
   use block_module, only: block_type
   use FD_module, only: apply_FDstencil
   implicit none

   private

   public :: full_weighing_restriction_2D, bilinear_prolongation_2D, apply_correction_2D

contains

   !> Full weighing restriction operator for 2D problems
   subroutine full_weighing_restriction_2D(fine_extended_dims, coarse_extended_dims, fine_residual, coarse_residual)
      integer, dimension(:), intent(in) :: fine_extended_dims, coarse_extended_dims
      real, contiguous, dimension(:,:), intent(in) :: fine_residual
      real, dimension(:,:), intent(out) :: coarse_residual

      real, dimension(3,3), parameter :: full_weighing_restriction_2D_stencil = &
         reshape([1.0, 2.0, 1.0, &
         2.0, 4.0, 2.0, &
         1.0, 2.0, 1.0]/16.0, [3,3])

      integer :: ii, jj
      integer :: fi, fj

      ! Initialize the coarse solution with zeros if the entire array is not guaranteed to be filled. This step might not be necessary...

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(fine_extended_dims, coarse_extended_dims, fine_residual, coarse_residual) &
      !$omp private(ii, jj, fi, fj)
      do ii = 1, coarse_extended_dims(2)
         do jj = 1, coarse_extended_dims(1)
            coarse_residual(jj,ii) = 0.0
         end do
      end do
      !$omp end parallel do

      ! Restrict the residual to the coarse grid

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(fine_extended_dims, coarse_extended_dims, fine_residual, coarse_residual) &
      !$omp private(ii, jj, fi, fj)
      do ii = 1, coarse_extended_dims(2)
         do jj = 1, coarse_extended_dims(1)
            fi = min(2 * ii, fine_extended_dims(2)-1)
            fj = min(2 * jj, fine_extended_dims(1)-1)

            coarse_residual(jj,ii) = sum(fine_residual(fj-1:fj+1, fi-1:fi+1) * full_weighing_restriction_2D_stencil)

         end do
      end do
      !$omp end parallel do

   end subroutine full_weighing_restriction_2D

   !> Bilinear prolongation operator for 2D problems
   subroutine bilinear_prolongation_2D(fine_extended_dims, coarse_extended_dims, fine_residual, coarse_solution)
      integer, dimension(:), intent(in) :: fine_extended_dims, coarse_extended_dims
      real, dimension(:,:), intent(out) :: fine_residual
      real, dimension(:,:), intent(in) :: coarse_solution

      integer :: ii, jj
      integer :: ci, cj
      real :: f00, f01, f10, f11

      ! Initialize the fine solution with zeros if the entire array is not guaranteed to be filled.

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(fine_extended_dims, fine_residual) &
      !$omp private(ii, jj)
      do ii = 1, fine_extended_dims(2)
         do jj = 1, fine_extended_dims(1)
            fine_residual(jj,ii) = 0.0
         end do
      end do
      !$omp end parallel do

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(coarse_extended_dims, fine_extended_dims, coarse_solution, fine_residual) &
      !$omp private(ii, jj, ci, cj, f00, f01, f10, f11)
      do ii = 1, coarse_extended_dims(2)-1
         do jj = 1, coarse_extended_dims(1)-1
            ci = min(2*ii, fine_extended_dims(2)-1)
            cj = min(2*jj, fine_extended_dims(1)-1)

            f00 = coarse_solution(jj, ii) ! f(j,i)
            f10 = coarse_solution(jj+1, ii) ! f(j+1,i)
            f01 = coarse_solution(jj, ii+1) ! f(j,i+1)
            f11 = coarse_solution(jj+1, ii+1) ! f(j+1,i+1)

            fine_residual(cj, ci) = (f00 + f10 + f01 + f11) / 4.0
            fine_residual(cj-1, ci) = (f00 + f10) / 2.0
            fine_residual(cj, ci-1) = (f00 + f01) / 2.0
            fine_residual(cj-1, ci-1) = f00

            if(ci == fine_extended_dims(2)-1) then
               fine_residual(cj, ci+1) = (f01 + f11) / 2.0
               fine_residual(cj-1, ci+1) = f01
            end if

            if(cj == fine_extended_dims(1)-1) then
               fine_residual(cj+1, ci) = (f10 + f11) / 2.0
               fine_residual(cj+1, ci-1) = f10
            end if

         end do
      end do
      !$omp end parallel do

   end subroutine bilinear_prolongation_2D

   !> Subroutine to apply the error correction to the solution
   subroutine apply_correction_2D(fine_extended_dims, fine_solution, correction)
      integer, dimension(:), intent(in) :: fine_extended_dims
      real, dimension(:,:), intent(inout) :: fine_solution
      real, dimension(:,:), intent(in) :: correction ! Stored in the residual array

      integer :: fi, fj

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(fine_extended_dims, fine_solution, correction) &
      !$omp private(fi, fj)
      do fi = 1, fine_extended_dims(2)
         do fj = 1, fine_extended_dims(1)
            fine_solution(fj,fi) = fine_solution(fj,fi) + correction(fj,fi)
         end do
      end do
      !$omp end parallel do

   end subroutine apply_correction_2D

end module multigrid_module
