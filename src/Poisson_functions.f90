module Poisson_functions_module
   use utility_functions_module, only: reshape_real_1D_to_2D, set_bc_zero_2D
   use block_module, only: block_type
   use FD_module, only: FDstencil_type, calculate_scaled_coefficients, get_coefficients_wrapper, &
      apply_FDstencil_2D, update_value_from_stencil_2D, set_matrix_coefficients
   implicit none

   private

   public :: Poisson_Gauss_Seidel_2D, Poisson_Gauss_Seidel_RB_2D, Poisson_Jacobi_2D, Poisson_assemble_matrix_2D
   public :: Poisson_residual_2D

contains

   !> The Gauss-Seidel iteration for the 2D Poisson problem fully sequential can use all stencils.
   !! Make sure the stencils has been scaled before calling this function
   subroutine Poisson_Gauss_Seidel_2D(omega, dxx_loc, dyy_loc, data_block, FDstencil, &
      local_u_diff_norm, local_r_diff_norm, local_u_norm, local_r_norm)
      real, intent(in) :: omega
      integer, intent(in) :: dxx_loc, dyy_loc
      type(block_type), intent(inout) :: data_block
      type(FDstencil_type), intent(inout) :: FDstencil
      real, intent(out) :: local_u_diff_norm, local_r_diff_norm, local_u_norm, local_r_norm

      integer :: ii, jj
      integer, dimension(2) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
      real :: f_val, u0_val, u1_val, r0_val, r1_val

      local_u_diff_norm = 0.0
      local_r_diff_norm = 0.0
      local_u_norm = 0.0
      local_r_norm = 0.0

      ! Sequential loop
      do jj = data_block%block_begin_c(2)+1, data_block%block_end_c(2)
         do ii = data_block%block_begin_c(1)+1, data_block%block_end_c(1)
            local_indices = [ii,jj]

            f_val = data_block%f_matrix_ptr_2D(local_indices(1),local_indices(2))
            u0_val = data_block%matrix_ptr_2D(local_indices(1),local_indices(2))

            call get_coefficients_wrapper(FDstencil, [1,1], data_block%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(dxx_loc * FDstencil%num_stencil_elements + 1:(dxx_loc + 1) * FDstencil%num_stencil_elements)
            dyy => coefficients(dyy_loc * FDstencil%num_stencil_elements + 1:(dyy_loc + 1) * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dyy

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

            call update_value_from_stencil_2D(combined_stencils_2D, data_block%matrix_ptr_2D, local_indices, alpha, beta, &
               f_val, u1_val, r1_val)

            u1_val = (1.0 - omega) * u0_val + omega * u1_val

            data_block%matrix_ptr_2D(local_indices(1),local_indices(2)) = u1_val
            local_u_diff_norm = local_u_diff_norm + (abs(u1_val - u0_val)**2)
            local_u_norm = local_u_norm + (abs(u1_val)**2)

         end do
      end do

   end subroutine Poisson_Gauss_Seidel_2D

   !> The Red-Black Gauss-Seidel iteration for the 2D Poisson problem but only for the 5-point stencil.
   subroutine Poisson_Gauss_Seidel_RB_2D(omega, dxx_loc, dyy_loc, data_block, FDstencil, &
      local_u_diff_norm, local_r_diff_norm, local_u_norm, local_r_norm)
      real, intent(in) :: omega
      integer, intent(in) :: dxx_loc, dyy_loc
      type(block_type), intent(inout) :: data_block
      type(FDstencil_type), intent(inout) :: FDstencil
      real, intent(out) :: local_u_diff_norm, local_r_diff_norm, local_u_norm, local_r_norm

      integer :: ii, jj, color
      integer, dimension(2) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
      real :: f_val, u0_val, u1_val, r0_val, r1_val

      local_u_diff_norm = 0.0
      local_r_diff_norm = 0.0
      local_u_norm = 0.0
      local_r_norm = 0.0

      ! Red-Black loop
      do color = 0, 1
         !$omp parallel do collapse(2) default(none) reduction(+:local_u_diff_norm, local_u_norm) &
         !$omp shared(omega, dxx_loc, dyy_loc, data_block, FDstencil, color) &
         !$omp private(ii, jj, local_indices, alpha, beta, coefficients, dxx, dyy, combined_stencils, &
         !$omp combined_stencils_2D, f_val, u0_val, u1_val, r0_val, r1_val)
         do jj = data_block%block_begin_c(2)+1, data_block%block_end_c(2)
            do ii = data_block%block_begin_c(1)+1, data_block%block_end_c(1)
               if(mod(ii+jj,2) /= color) cycle ! Avoid the if-statement somehow...
               local_indices = [ii,jj]

               f_val = data_block%f_matrix_ptr_2D(local_indices(1),local_indices(2))
               u0_val = data_block%matrix_ptr_2D(local_indices(1),local_indices(2))

               call get_coefficients_wrapper(FDstencil, [1,1], data_block%extended_block_dims, local_indices, &
                  alpha, beta, coefficients)

               dxx => coefficients(dxx_loc * FDstencil%num_stencil_elements + 1:(dxx_loc + 1) * FDstencil%num_stencil_elements)
               dyy => coefficients(dyy_loc * FDstencil%num_stencil_elements + 1:(dyy_loc + 1) * FDstencil%num_stencil_elements)

               combined_stencils = dxx + dyy

               call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

               call update_value_from_stencil_2D(combined_stencils_2D, data_block%matrix_ptr_2D, local_indices, alpha, beta, &
                  f_val, u1_val, r1_val)

               u1_val = (1.0 - omega) * u0_val + omega * u1_val

               data_block%matrix_ptr_2D(local_indices(1),local_indices(2)) = u1_val
               local_u_diff_norm = local_u_diff_norm + (abs(u1_val - u0_val)**2)
               local_u_norm = local_u_norm + (abs(u1_val)**2)

            end do
         end do

         !$omp end parallel do

      end do

   end subroutine Poisson_Gauss_Seidel_RB_2D

   !> The Jacobi iteration for the 2D Poisson problem fully parallelized for all stencils.
   !! Make sure the stencils has been scaled before calling this function
   subroutine Poisson_Jacobi_2D(omega, dxx_loc, dyy_loc, data_block, FDstencil, &
      local_u_diff_norm, local_r_diff_norm, local_u_norm, local_r_norm)
      real, intent(in) :: omega
      integer, intent(in) :: dxx_loc, dyy_loc
      type(block_type), intent(inout) :: data_block
      type(FDstencil_type), intent(inout) :: FDstencil
      real, intent(out) :: local_u_diff_norm, local_r_diff_norm, local_u_norm, local_r_norm

      integer :: ii, jj
      integer, dimension(2) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
      real :: f_val, u0_val, u1_val, r0_val, r1_val

      local_u_diff_norm = 0.0
      local_r_diff_norm = 0.0
      local_u_norm = 0.0
      local_r_norm = 0.0

      !$omp parallel do collapse(2) default(none) reduction(+:local_u_diff_norm, local_u_norm) &
      !$omp shared(omega, dxx_loc, dyy_loc, data_block, FDstencil) &
      !$omp private(ii, jj, local_indices, alpha, beta, coefficients, dxx, dyy, combined_stencils, &
      !$omp combined_stencils_2D, f_val, u0_val, u1_val, r0_val, r1_val)
      do jj = data_block%block_begin_c(2)+1, data_block%block_end_c(2)
         do ii = data_block%block_begin_c(1)+1, data_block%block_end_c(1)
            local_indices = [ii,jj]

            f_val = data_block%f_matrix_ptr_2D(local_indices(1),local_indices(2))
            u0_val = data_block%matrix_ptr_2D(local_indices(1),local_indices(2))

            call get_coefficients_wrapper(FDstencil, [1,1], data_block%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(dxx_loc * FDstencil%num_stencil_elements + 1:(dxx_loc + 1) * FDstencil%num_stencil_elements)
            dyy => coefficients(dyy_loc * FDstencil%num_stencil_elements + 1:(dyy_loc + 1) * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dyy

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

            call update_value_from_stencil_2D(combined_stencils_2D, data_block%matrix_ptr_2D, local_indices, alpha, beta, &
               f_val, u1_val, r1_val)

            u1_val = (1.0 - omega) * u0_val + omega * u1_val

            data_block%temp_matrix_ptr_2D(local_indices(1),local_indices(2)) = u1_val
            local_u_diff_norm = local_u_diff_norm + (abs(u1_val - u0_val)**2)
            local_u_norm = local_u_norm + (abs(u1_val)**2)

         end do
      end do

      !$omp end parallel do

   end subroutine Poisson_Jacobi_2D

   !> The global assembly of the matrix A for the 2D Poisson problem fully parallelized for all stencils.
   subroutine Poisson_assemble_matrix_2D(dxx_loc, dyy_loc, data_block, FDstencil)
      integer, intent(in) :: dxx_loc, dyy_loc
      type(block_type), intent(inout) :: data_block
      type(FDstencil_type), intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(2) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils

      !> Calculate the scaled coefficients for the stencils
      call calculate_scaled_coefficients(data_block%ndims, data_block%extended_grid_dx, FDstencil)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(dxx_loc, dyy_loc, data_block, FDstencil) &
      !$omp private(ii, jj, local_indices, alpha, beta, coefficients, dxx, dyy, combined_stencils)
      do jj = data_block%extended_block_begin_c(2)+1, data_block%extended_block_end_c(2)
         do ii = data_block%extended_block_begin_c(1)+1, data_block%extended_block_end_c(1)
            local_indices = [ii,jj]

            call get_coefficients_wrapper(FDstencil, [1,1], data_block%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(dxx_loc * FDstencil%num_stencil_elements + 1:(dxx_loc + 1) * FDstencil%num_stencil_elements)
            dyy => coefficients(dyy_loc * FDstencil%num_stencil_elements + 1:(dyy_loc + 1) * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dyy

            call set_matrix_coefficients(data_block%ndims, FDstencil%stencil_sizes, data_block%extended_block_dims, &
               combined_stencils, data_block%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)

         end do
      end do

      !$omp end parallel do

   end subroutine Poisson_assemble_matrix_2D

   !> Calculate the residual of the 2D-Poisson equation with zero boundaries.
   subroutine Poisson_residual_2D(dxx_loc, dyy_loc, data_block, FDstencil)
      integer, intent(in) :: dxx_loc, dyy_loc
      type(block_type), intent(inout) :: data_block
      type(FDstencil_type), intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(2) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
      real :: laplacian_p

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(dxx_loc, dyy_loc, data_block, FDstencil) &
      !$omp private(ii, jj, local_indices, alpha, beta, coefficients, dxx, dyy, combined_stencils, &
      !$omp combined_stencils_2D, laplacian_p)
      do jj = data_block%block_begin_c(2)+1, data_block%block_end_c(2)
         do ii = data_block%block_begin_c(1)+1, data_block%block_end_c(1)
            local_indices = [ii,jj]

            call get_coefficients_wrapper(FDstencil, [1,1], data_block%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(dxx_loc * FDstencil%num_stencil_elements + 1:(dxx_loc + 1) * FDstencil%num_stencil_elements)
            dyy => coefficients(dyy_loc * FDstencil%num_stencil_elements + 1:(dyy_loc + 1) * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dyy

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

            call apply_FDstencil_2D(combined_stencils_2D, data_block%matrix_ptr_2D, local_indices, alpha, beta, laplacian_p)

            data_block%residual_matrix_ptr_2D(local_indices(1),local_indices(2)) = &
               data_block%f_matrix_ptr_2D(local_indices(1),local_indices(2)) - laplacian_p

         end do
      end do

      !$omp end parallel do

      ! Depends on the boundary conditions this is for Dirichlet. For Neumann we need to calculate the gradient at the point.
      call set_bc_zero_2D(data_block%residual_matrix_ptr_2D)

   end subroutine Poisson_residual_2D

end module Poisson_functions_module
