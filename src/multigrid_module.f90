!> This module contains implementations of the necessary functions for the multigrid method
module multigrid_module
   use utility_functions_module, only: IDX_XD, IDX_XD_INV
   use block_module, only: block_type
   use FD_module, only: apply_FDstencil
   implicit none

   private

   public :: full_weighing_restriction_2D, nearest_neighbor_prolongation_2D

contains

   !> Full weighing restriction operator for 2D problems
   subroutine full_weighing_restriction_2D(fine_extended_dims, coarse_extended_dims, fine_residual, coarse_residual)
      integer, dimension(:), intent(in) :: fine_extended_dims, coarse_extended_dims
      real, dimension(:,:), intent(in) :: fine_residual
      real, dimension(:,:), intent(out) :: coarse_residual

      integer :: i, j

      ! Make it parallel
      do i = 1, coarse_extended_dims(1)
         do j = 1, coarse_extended_dims(2)
            coarse_residual(i, j) = (4.0 * fine_residual(2*i, 2*j) + &
               2.0 * (fine_residual(2*i-1, 2*j) + fine_residual(2*i+1, 2*j) + &
               fine_residual(2*i, 2*j-1) + fine_residual(2*i, 2*j+1)) + &
               1.0 * (fine_residual(2*i-1, 2*j-1) + fine_residual(2*i+1, 2*j-1) + &
               fine_residual(2*i-1, 2*j+1) + fine_residual(2*i+1, 2*j+1))) / 16.0
         end do
      end do

   end subroutine full_weighing_restriction_2D

   !> Nearest neighbor prolongation operator for 2D problems
   subroutine nearest_neighbor_prolongation_2D(coarse_extended_dims, fine_extended_dims, coarse_solution, fine_solution)
      integer, dimension(:), intent(in) :: coarse_extended_dims, fine_extended_dims
      real, dimension(:,:), intent(in) :: coarse_solution
      real, dimension(:,:), intent(out) :: fine_solution

      integer :: i, j

      ! Initialize the fine solution with zeros
      fine_solution = 0.0

      ! Perform nearest neighbor prolongation
      do i = 1, coarse_extended_dims(1)
         do j = 1, coarse_extended_dims(2)
            fine_solution(2*i-1, 2*j-1) = coarse_solution(i, j)
            fine_solution(2*i-1, 2*j) = coarse_solution(i, j)
            fine_solution(2*i, 2*j-1) = coarse_solution(i, j)
            fine_solution(2*i, 2*j) = coarse_solution(i, j)
         end do
      end do

   end subroutine nearest_neighbor_prolongation_2D

end module multigrid_module
