module TravelingWave_Poisson_2D_module
   use constants_module, only: pi, MASTER_RANK
   use mpi, only: MPI_WTIME, MPI_DOUBLE_PRECISION, MPI_SUM
   use utility_functions_module, only:open_txt_file, close_txt_file, sleeper_function, &
      reshape_real_1D_to_2D
   use mpi_wrapper_module, only: write_block_data_to_file, all_reduce_mpi_wrapper
   use solver_module, only: SolverParamsType, set_SolverParamsType, ResultType, print_resultType
   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, print_cart_comm_type
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type, &
      sendrecv_data_neighbors
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      calculate_scaled_coefficients, get_coefficients_wrapper, &
      update_value_from_stencil_2D, apply_FDstencil_2D, apply_FDstencil, update_value_from_stencil
   use functions_module, only: calculate_point
   use solver_module, only: check_convergence
   use multigrid_module, only: full_weighing_restriction_2D, bilinear_prolongation_2D
   implicit none

   private

   !> Number of dimensions and number of derivatives
   integer, parameter :: ndims = 2, num_derivatives = 4
   integer, dimension(ndims * num_derivatives), parameter :: derivatives = [1,0,0,1,2,0,0,2] ! ds, dx, dss, dxx

   !> Physical parameters
   real, parameter :: g = 9.81, mu = 1.3059*(1e-6), rho = 1000.0, nu = mu/rho

   !> Domain parameters
   real, parameter :: Ls = 1.0, Lx = 1.0

   !> Wave parameters
   real, parameter :: start_t = 0.0
   real, parameter :: kh = 1.0, lwave = Lx, kwave = 2.0*pi/lwave, hd = kh/kwave, cwave = sqrt(g/kwave*tanh(kh)), &
      Twave = lwave/cwave, wwave=2.0*pi/Twave, Hwave = 0.02

   !> Grid parameters
   integer, dimension(ndims), parameter :: grid_size = [64,64], processor_dims = [1,1]
   logical, dimension(ndims), parameter :: periods = [.false., .true.]
   logical, parameter :: reorder = .true.
   real, dimension(ndims), parameter :: domain_begin = [0,0], domain_end = [Ls,Lx]
   integer, dimension(ndims), parameter :: stencil_sizes = 3
   integer, dimension(ndims), parameter :: uv_ghost_begin = [0,0], uv_ghost_end = [0,0]
   integer, dimension(ndims), parameter :: p_ghost_begin = [1,0], p_ghost_end = [1,0]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   !> Solver parameters
   real, parameter :: tol = 1e-16, div_tol = 1e-1, omega = 1.0
   integer, parameter :: max_iter = 10000*(1.0 + 1.0 - omega), multigrid_max_level = 1

   public :: TravelingWave_Poisson_2D_main

contains

   subroutine TravelingWave_Poisson_2D_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      type(SolverParamsType) :: solver_params
      type(ResultType) :: result_params
      type(comm_type) :: comm_params
      type(block_type) :: u_params, v_params, p_params
      type(FDstencil_type) :: FDstencil_params

      real, dimension(4) :: result_array_with_timings
      integer :: iounit

      !call sleeper_function(1)

      ! Set the solver parameters: tol, div_tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(tol, div_tol, max_iter, 2, solver_params)

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm_params)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, u_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, v_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         p_ghost_begin, p_ghost_end, stencil_begin, stencil_end, p_params)

      call Write_analytical_function(u_params, v_params, p_params)

      call sendrecv_data_neighbors(comm_params%comm, u_params, u_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, v_params, v_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, p_params, p_params%matrix_ptr)

      ! Time and solve the system
      result_array_with_timings(1) = MPI_WTIME()

      ! Write the velocity convergence to the f-matrix of the pressure block
      call calculate_velocity_divergence(comm_params, u_params, v_params, FDstencil_params, p_params)

      ! Run a W-cycle
      call w_cycle(solver_params, comm_params, p_params, FDstencil_params, result_params, multigrid_max_level, 1)

      call calculate_residual(comm_params, p_params, FDstencil_params)

      result_array_with_timings(2) = MPI_WTIME()
      result_array_with_timings(3) = result_array_with_timings(2) - result_array_with_timings(1)
      call all_reduce_mpi_wrapper(result_array_with_timings(3), result_array_with_timings(4), 1, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_params%comm)

      ! Write out the result and timings from the master rank
      if(rank == MASTER_RANK) then
         write(*,*) "Results: "
         call print_resultType(result_params)
         write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(4)/world_size, " seconds"
      end if

      call sendrecv_data_neighbors(comm_params%comm, u_params, u_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, v_params, v_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, p_params, p_params%matrix_ptr)

      ! Write the block data to a .txt file
      call open_txt_file("output/output_from_", rank, iounit)
      call print_cart_comm_type(comm_params, iounit)
      !call print_block_type(u_params, iounit)
      !call print_block_type(v_params, iounit)
      call print_block_type(p_params, iounit)
      call close_txt_file(iounit)

      ! Write the block data to a .dat file
      call write_block_data_to_file(u_params%data_layout, "output/u_solution.dat", comm_params%comm, u_params%matrix_ptr)
      call write_block_data_to_file(v_params%data_layout, "output/v_solution.dat", comm_params%comm, v_params%matrix_ptr)
      call write_block_data_to_file(p_params%data_layout, "output/p_solution.dat", comm_params%comm, p_params%matrix_ptr)
      call write_block_data_to_file(p_params%data_layout, "output/vel_div.dat", comm_params%comm, p_params%f_matrix_ptr)
      call write_block_data_to_file(p_params%data_layout, "output/residual.dat", comm_params%comm, p_params%residual_matrix_ptr)

      ! Deallocate used memory

      call deallocate_cart_comm_type(comm_params)

      call deallocate_finite_difference_stencil(FDstencil_params)

      call deallocate_block_type(u_params)

      call deallocate_block_type(v_params)

      call deallocate_block_type(p_params)

   end subroutine TravelingWave_Poisson_2D_main

   recursive subroutine w_cycle(solver, comm, p_block_fine, FDstencil, result, max_level, level)
      type(SolverParamsType), intent(in) :: solver
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: p_block_fine
      type(FDstencil_type), intent(inout) :: FDstencil
      type(ResultType), intent(out) :: result
      integer, intent(in) :: max_level, level

      type(block_type) :: p_block_coarse

      write(0,*) "Multigrid level: ", level

      ! Pre-smoothing
      call solve_pressure_poisson(comm, p_block_fine, FDstencil, solver, result)

      ! If max level is not reached we restrict the residual to a coarser grid
      if(level /= max_level) then

         ! Create a coarser grid
         call create_block_type(ndims, 1, 1, domain_begin, domain_end, p_block_fine%grid_size/2, comm, &
            p_ghost_begin, p_ghost_end, stencil_begin, stencil_end, p_block_coarse)

         ! Compute residual errors
         call calculate_residual(comm, p_block_fine, FDstencil)

         ! Restrict the residual to the coarser grid
         call full_weighing_restriction_2D(p_block_fine%extended_block_dims, p_block_coarse%extended_block_dims, &
            p_block_fine%residual_matrix_ptr_2D, p_block_coarse%f_matrix_ptr_2D)

         ! We cycle since we are not at the max level
         call w_cycle(solver, comm, p_block_coarse, FDstencil, result, max_level, level+1)

         ! Prolongate the solution to the finer grid
         call bilinear_prolongation_2D(p_block_fine%extended_block_dims, p_block_coarse%extended_block_dims, &
            p_block_fine%residual_matrix_ptr_2D, p_block_coarse%matrix_ptr_2D)

         ! Error correction
         p_block_fine%matrix_ptr_2D = p_block_fine%matrix_ptr_2D + p_block_fine%residual_matrix_ptr_2D

         ! Second-smoothing
         call solve_pressure_poisson(comm, p_block_fine, FDstencil, solver, result)

         ! Calculate the residual again
         call calculate_residual(comm, p_block_fine, FDstencil)

         ! Restrict the residual to the coarser grid again
         call full_weighing_restriction_2D(p_block_fine%extended_block_dims, p_block_coarse%extended_block_dims, &
            p_block_fine%residual_matrix_ptr_2D, p_block_coarse%f_matrix_ptr_2D)

         ! Second cycle
         call w_cycle(solver, comm, p_block_coarse, FDstencil, result, max_level, level+1)

         ! Prolongate the solution to the finer grid
         call bilinear_prolongation_2D(p_block_fine%extended_block_dims, p_block_coarse%extended_block_dims, &
            p_block_fine%residual_matrix_ptr_2D, p_block_coarse%matrix_ptr_2D)

         ! Final error correction
         p_block_fine%matrix_ptr_2D = p_block_fine%matrix_ptr_2D + p_block_fine%residual_matrix_ptr_2D

         ! Post-smoothing
         call solve_pressure_poisson(comm, p_block_fine, FDstencil, solver, result)

         ! Deallocate the coarser grid
         call deallocate_block_type(p_block_coarse)

      end if

   end subroutine w_cycle

   ! Calculate the residual of the pressure poisson equation
   subroutine calculate_residual(comm, p_block, FDstencil)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: p_block
      type(FDstencil_type), intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(p_block%ndims) :: p_local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dss
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
      real :: laplacian_p

      call calculate_scaled_coefficients(p_block%ndims, p_block%extended_grid_dx, FDstencil)

      call sendrecv_data_neighbors(comm%comm, p_block, p_block%matrix_ptr)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(p_block, FDstencil) &
      !$omp private(ii, jj, p_local_indices, alpha, beta, coefficients, dxx, dss, combined_stencils, &
      !$omp combined_stencils_2D, laplacian_p)
      do ii = p_block%block_begin_c(2)+1, p_block%block_end_c(2)
         do jj = p_block%block_begin_c(1)+1, p_block%block_end_c(1)
            p_local_indices = [jj,ii]

            call get_coefficients_wrapper(FDstencil, p_block%block_begin_c+1, p_block%block_end_c, &
               p_local_indices, alpha, beta, coefficients)

            dss => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
            dxx => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dss

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)
            call apply_FDstencil_2D(combined_stencils_2D, p_block%matrix_ptr_2D, p_local_indices, alpha, beta, laplacian_p)

            p_block%residual_matrix_ptr_2D(jj,ii) = p_block%f_matrix_ptr_2D(jj,ii) - laplacian_p

         end do
      end do

      !$omp end parallel do

   end subroutine calculate_residual

   ! Calculate velocity divergence
   subroutine calculate_velocity_divergence(comm, u_block, v_block, FDstencil, p_block)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(u_block%ndims) :: uv_local_indices, p_local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dx, ds
      real, contiguous, dimension(:,:), pointer :: dx_2D, ds_2D
      real :: u_x, v_s

      call calculate_scaled_coefficients(u_block%ndims, u_block%extended_grid_dx, FDstencil)

      call sendrecv_data_neighbors(comm%comm, u_block, u_block%matrix_ptr)
      call sendrecv_data_neighbors(comm%comm, v_block, v_block%matrix_ptr)
      call sendrecv_data_neighbors(comm%comm, p_block, p_block%matrix_ptr)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(u_block, v_block, FDstencil, p_block) &
      !$omp private(ii, jj, uv_local_indices, p_local_indices, alpha, beta, coefficients, dx, ds, dx_2D, ds_2D, u_x, v_s)
      do ii = u_block%block_begin_c(2)+1, u_block%block_end_c(2)
         do jj = u_block%block_begin_c(1)+1, u_block%block_end_c(1)
            uv_local_indices = [jj,ii]
            p_local_indices = p_block%block_begin_c + uv_local_indices - u_block%block_begin_c

            call get_coefficients_wrapper(FDstencil, u_block%block_begin_c+1, u_block%block_end_c, &
               uv_local_indices, alpha, beta, coefficients)

            ds => coefficients(1:FDstencil%num_stencil_elements)
            dx => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, ds, ds_2D)

            call apply_FDstencil_2D(dx_2D, u_block%matrix_ptr_2D, uv_local_indices, alpha, beta, u_x)
            call apply_FDstencil_2D(ds_2D, v_block%matrix_ptr_2D, uv_local_indices, alpha, beta, v_s)

            p_block%f_matrix_ptr_2D(p_local_indices(1), p_local_indices(2)) = (u_x + v_s)

         end do
      end do

      !$omp end parallel do

   end subroutine calculate_velocity_divergence

   subroutine solve_pressure_poisson(comm, p_block, FDstencil, solver, result)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: p_block
      type(FDstencil_type), intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver
      type(ResultType), intent(out) :: result

      integer :: it, converged
      real, dimension(4) :: norm_array

      call calculate_scaled_coefficients(p_block%ndims, p_block%extended_grid_dx, FDstencil)

      call write_pressure_BC(p_block%extended_block_dims, p_block%matrix_ptr_2D)
      call sendrecv_data_neighbors(comm%comm, p_block, p_block%matrix_ptr)

      norm_array = [1e3,1e6,1e9,1e12]

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < solver%max_iter)
         call poisson_iteration(p_block, FDstencil, norm_array(1))

         call write_pressure_BC(p_block%extended_block_dims, p_block%matrix_ptr_2D)
         call sendrecv_data_neighbors(comm%comm, p_block, p_block%matrix_ptr)

         call check_convergence(comm%comm, solver%tol, solver%divergence_tol, &
            1.0/product(p_block%extended_grid_size), norm_array, converged)
         result%converged = converged
         if(converged == -1 .and. it > 0) then
            write(*,*) "Convergence failed"
            exit
         end if

         it = it + 1

      end do

      result%iterations = it
      result%global_norm = norm_array(3)
      result%relative_norm = norm_array(4)

   end subroutine solve_pressure_poisson

   subroutine poisson_iteration(p_block, FDstencil, local_norm)
      type(block_type), intent(inout) :: p_block
      type(FDstencil_type), intent(inout) :: FDstencil
      real, intent(out) :: local_norm

      integer :: color, ii, jj
      integer, dimension(ndims) :: p_local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dss
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
      real :: old_val, f_val, new_val

      local_norm = 0.0

      ! Red-black Gauss-Seidel
      do color = 0, 1
         !$omp parallel do collapse(2) default(none) reduction(+:local_norm) &
         !$omp shared(p_block, FDstencil, color) &
         !$omp private(ii, jj, p_local_indices, alpha, beta, coefficients, dxx, dss, combined_stencils, &
         !$omp combined_stencils_2D, old_val, f_val, new_val)
         do ii = p_block%block_begin_c(2)+1, p_block%block_end_c(2)
            do jj = p_block%block_begin_c(1)+1, p_block%block_end_c(1)
               if(mod(ii+jj,2) /= color) cycle ! Avoid the if-statement somehow...
               p_local_indices = [jj,ii]

               old_val = p_block%matrix_ptr_2D(jj,ii)
               f_val = p_block%f_matrix_ptr_2D(jj,ii)

               call get_coefficients_wrapper(FDstencil, [1,1], p_block%extended_block_dims, &
                  p_local_indices, alpha, beta, coefficients)

               dss => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
               dxx => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

               combined_stencils = dxx + dss

               call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

               call update_value_from_stencil_2D(combined_stencils_2D, p_block%matrix_ptr_2D, p_local_indices, alpha, beta, &
                  f_val, new_val)

               new_val = (1.0 - omega) * old_val + omega * new_val

               p_block%matrix_ptr_2D(jj,ii) = new_val

               local_norm = local_norm + (abs(new_val - old_val)**2)

            end do
         end do

         !$omp end parallel do

      end do

   end subroutine poisson_iteration

   !> Write the boundary conditions for the pressure. Pressure is set to zero at the top and P_y = 0 at the bottom.
   !! Periodic boundary conditions are used for the left and right walls.
   !! Currently this is wrong because we use the local extended dims. But since we send and recv the data is overwritten everywhere where it is wrong. Should fix it though!
   subroutine write_pressure_BC(dims, matrix)
      integer, dimension(:), intent(in) :: dims
      real, dimension(:,:), intent(inout) :: matrix

      integer :: ii

      ! Periodic is built into the matrix, so we don't need to do anything here.

      !$omp parallel do default(none) &
      !$omp shared(dims, matrix) &
      !$omp private(ii)
      do ii = 1, dims(2)
         matrix(1,ii) = matrix(3,ii) ! At seabottom: p_y = 0 (Neumann condition)
         matrix(dims(1),ii) = matrix(dims(1)-2,ii)! At the free surface: p = 0 (Dirichlet condition) ! Replacing with a neumann condition improves the error in this part! matrix(dims(1),dims(2)) = 0.0
      end do

      !$omp end parallel do

   end subroutine write_pressure_BC

   subroutine Write_analytical_function(u_block, v_block, p_block)
      type(block_type), intent(inout) :: u_block, v_block, p_block

      integer :: ii, jj
      integer, dimension(u_block%ndims) :: uv_local_indices, p_local_indices, uv_global_indices
      real, dimension(u_block%ndims) :: point
      real :: x, sigma
      real :: uu, uu_x, ww, ww_z, eta, etax, etaxx, pp

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(u_block, v_block, p_block) &
      !$omp private(ii, jj, uv_local_indices, p_local_indices, uv_global_indices, &
      !$omp point, x, sigma, uu, uu_x, ww, ww_z, eta, etax, etaxx, pp)
      do ii = u_block%extended_global_begin_c(2)+1, u_block%extended_block_end_c(2)
         do jj = u_block%extended_global_begin_c(1)+1, u_block%extended_block_end_c(1)
            uv_local_indices = [jj,ii]
            p_local_indices = p_block%block_begin_c + uv_local_indices - u_block%block_begin_c
            uv_global_indices = u_block%extended_global_begin_c + uv_local_indices

            call calculate_point(u_block%ndims, -u_block%global_grid_begin, uv_global_indices, &
               u_block%domain_begin, u_block%grid_dx, point)

            sigma = point(1) ! we should do 1-y to start from the surface! Also do it in python code
            x = point(2)

            call TravelingWave1D(Hwave, cwave, kwave, sigma, hd, wwave, start_t, x, uu, uu_x, ww, ww_z, eta, etax, etaxx, pp)

            u_block%matrix_ptr_2D(uv_local_indices(1), uv_local_indices(2)) = uu
            v_block%matrix_ptr_2D(uv_local_indices(1), uv_local_indices(2)) = ww
            p_block%matrix_ptr_2D(p_local_indices(1), p_local_indices(2)) = 0.0!pp
            p_block%f_matrix_ptr_2D(p_local_indices(1), p_local_indices(2)) = 0.0!uu_x + ww_z

         end do
      end do

      !$omp end parallel do

   end subroutine Write_analytical_function

   elemental subroutine TravelingWave1D(H, c, k, sigma, h_small, w, t, x, uu, uu_x, ww, ww_z, eta, etax, etaxx, pp)
      real, intent(in)  :: H, c, k, sigma, h_small, w, t, x
      real, intent(out) :: uu, uu_x, ww, ww_z, eta, etax, etaxx, pp

      real :: z

      z = sigma*h_small - h_small

      ! Compute analytical solution at level z
      uu    =  k*H*c/2*cosh(k*(z + h_small))/sinh(k * h_small)*cos(w*t-k*x)
      uu_x  = k**2 * H * c / 2 * cosh(k * (z + h_small)) / sinh(k * h_small) * sin(w * t - k * x)

      ww    = -k*H*c/2*sinh(k*(z + h_small))/sinh(k*h_small)*sin(w*t-k*x)
      ww_z  = -k**2 * H * c / 2 * cosh(k * (z + h_small)) / sinh(k * h_small) * sin(w * t - k * x)

      eta   =  H/2*cos(w*t-k*x)   ! z = 0
      etax  =  k*H/2*sin(w*t-k*x)
      etaxx = -(k**2)*H/2*cos(w*t-k*x);

      pp    = -H*c/2*cosh(k*(z+h_small))/sinh(k*h_small)*sin(w*t-k*x); ! see p. 83 in Svendsen & Jonsson (2001)
   end subroutine TravelingWave1D

end module TravelingWave_Poisson_2D_module
