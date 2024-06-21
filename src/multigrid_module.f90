!> This module contains implementations of the necessary functions for the multigrid method
module multigrid_module
   use utility_functions_module, only: IDX_XD, IDX_XD_INV
   use block_module, only: block_type
   use FD_module, only: apply_FDstencil
   implicit none

   private

   public :: full_weighing_restriction_2D, bilinear_prolongation_2D

contains

   !> Full weighing restriction operator for 2D problems
   subroutine full_weighing_restriction_2D(fine_extended_dims, coarse_extended_dims, fine_residual, coarse_residual)
      integer, dimension(:), intent(in) :: fine_extended_dims, coarse_extended_dims
      real, dimension(:,:), intent(in) :: fine_residual
      real, dimension(:,:), intent(out) :: coarse_residual

      real, dimension(3,3), parameter :: full_weighing_restriction_2D_stencil = &
         reshape([1.0, 2.0, 1.0, &
         2.0, 4.0, 2.0, &
         1.0, 2.0, 1.0]/16.0, [3,3])

      integer :: ii, jj
      integer :: fi, fj

      ! Initialize the coarse solution with zeros if the entire array is not guaranteed to be filled. Maybe parallelize this?
      coarse_residual = 0.0

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(fine_extended_dims, coarse_extended_dims, fine_residual, coarse_residual) &
      !$omp private(ii, jj, fi, fj)
      do ii = 1, coarse_extended_dims(2)
         do jj = 1, coarse_extended_dims(1)
            ! Compute corresponding fine grid indices
            fi = min(2 * ii, fine_extended_dims(2)-1)
            fj = min(2 * jj, fine_extended_dims(1)-1)

            coarse_residual(jj,ii) = sum(fine_residual(fj-1:fj+1, fi-1:fi+1) * full_weighing_restriction_2D_stencil)

         end do
      end do

      !$omp end parallel do

   end subroutine full_weighing_restriction_2D

   !> Bilinear prolongation operator for 2D problems
   subroutine bilinear_prolongation_2D(fine_extended_dims, coarse_extended_dims, fine_solution, coarse_solution)
      integer, dimension(:), intent(in) :: fine_extended_dims, coarse_extended_dims
      real, dimension(:,:), intent(out) :: fine_solution
      real, dimension(:,:), intent(in) :: coarse_solution

      integer :: ii, jj
      integer :: fi, fj

      ! Initialize the fine solution with zeros if the entire array is not guaranteed to be filled. Maybe parallelize this?
      fine_solution = 0.0

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(coarse_extended_dims, fine_extended_dims, coarse_solution, fine_solution) &
      !$omp private(ii, jj, fi, fj)
      do ii = 1, coarse_extended_dims(2)-1
         do jj = 1, coarse_extended_dims(1)-1
            ! Compute corresponding fine grid indices
            fi = min(2 * ii, fine_extended_dims(2)-1)
            fj = min(2 * jj, fine_extended_dims(1)-1)

            ! Directly assign values to corresponding points
            fine_solution(fj, fi) = coarse_solution(jj, ii)

            ! Interpolate horizontally
            fine_solution(fj, fi+1) = 0.5 * (coarse_solution(jj, ii) + coarse_solution(jj, ii+1))

            ! Interpolate vertically
            fine_solution(fj+1, fi) = 0.5 * (coarse_solution(jj, ii) + coarse_solution(jj+1, ii))

            ! Interpolate diagonally
            fine_solution(fj+1, fi+1) = 0.25 * (coarse_solution(jj, ii) + coarse_solution(jj, ii+1) + &
               coarse_solution(jj+1, ii) + coarse_solution(jj+1, ii+1))

         end do
      end do

      !$omp end parallel do

   end subroutine bilinear_prolongation_2D

end module multigrid_module
