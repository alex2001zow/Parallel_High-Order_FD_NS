module TravelingWave_Poisson_2D_module
   use constants_module, only: pi, MASTER_RANK
   use mpi, only: MPI_WTIME, MPI_DOUBLE_PRECISION, MPI_SUM
   use utility_functions_module, only: IDX_XD, IDX_XD_INV, open_txt_file, close_txt_file, sleeper_function, &
      reshape_real_1D_to_2D
   use mpi_wrapper_module, only: write_block_data_to_file, all_reduce_mpi_wrapper
   use solver_module, only: SolverParamsType, set_SolverParamsType
   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, print_cart_comm_type
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type, &
      sendrecv_data_neighbors
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      calculate_scaled_coefficients, get_FD_coefficients_from_index, update_value_from_stencil_2D, apply_stencil_2D
   use functions_module, only: calculate_point
   use solver_module, only: check_convergence
   use multigrid_module, only: full_weighing_restriction_2D, nearest_neighbor_prolongation_2D
   implicit none

   private


   integer, parameter :: ndims = 2, num_derivatives = 4
   integer, dimension(ndims*num_derivatives), parameter :: derivatives = [0,1,1,0,0,2,2,0] ! dx, ds, dxx, dss

   !> Physical parameters
   real, parameter :: g = 9.81, mu = 1.3059*(1e-6), rho = 1000.0, nu = mu/rho

   !> Domain parameters
   real, parameter :: Lx = 31.0, Ls = 1.0

   !> Wave parameters
   real, parameter :: start_t = 0.0
   real, parameter :: kh = 1.0, lwave = Lx, kwave = 2.0*pi/lwave, hd = kh/kwave, cwave = sqrt(g/kwave*tanh(kh)), &
      Twave = lwave/cwave, wwave=2.0*pi/Twave, Hwave = 0.02

   !> Grid parameters
   integer, dimension(2), parameter :: grid_size = [64,64], processor_dims = [1,1]
   logical, dimension(2), parameter :: periods = [.true., .false.]
   logical, parameter :: reorder = .true.
   real, dimension(2), parameter :: domain_begin = [0,0], domain_end = [Lx,Ls]
   integer, dimension(2), parameter :: stencil_sizes = [3,3]
   integer, dimension(2), parameter :: uv_ghost_begin = [0,0], uv_ghost_end = [0,0]
   integer, dimension(2), parameter :: p_ghost_begin = [0,1], p_ghost_end = [0,1]
   integer, dimension(2), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   !> Solver parameters
   real, parameter :: tol = 1e-9, div_tol = 1e-1
   integer, parameter :: max_iter = 20000

   public :: TravelingWave_Poisson_2D_main

contains

   subroutine TravelingWave_Poisson_2D_main(rank, world_size)
      integer, intent(in) :: rank, world_size


      type(SolverParamsType) :: solver_params
      type(comm_type) :: comm_params
      type(block_type) :: u_params, v_params, vel_div_params, p_params
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
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, vel_div_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         p_ghost_begin, p_ghost_end, stencil_begin, stencil_end, p_params)

      u_params%matrix_ptr = 0.0
      v_params%matrix_ptr = 0.0
      vel_div_params%matrix_ptr = 0.0
      p_params%matrix_ptr = 0.0

      call Write_analytical_function(u_params, v_params)

      call sendrecv_data_neighbors(comm_params%comm, u_params, u_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, v_params, v_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, p_params, p_params%matrix_ptr)

      ! Time and solve the system
      result_array_with_timings(1) = MPI_WTIME()

      call calculate_velocity_divergence(u_params, v_params, FDstencil_params, vel_div_params)

      call solve_pressure_poisson(comm_params, p_params, FDstencil_params, vel_div_params, solver_params)

      result_array_with_timings(2) = MPI_WTIME()
      result_array_with_timings(3) = result_array_with_timings(2) - result_array_with_timings(1)
      call all_reduce_mpi_wrapper(result_array_with_timings(3), result_array_with_timings(4), 1, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_params%comm)

      ! Write out the result and timings from the master rank
      if(rank == MASTER_RANK) then
         write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(4)/world_size, " seconds"
      end if

      call sendrecv_data_neighbors(comm_params%comm, u_params, u_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, v_params, v_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, p_params, p_params%matrix_ptr)

      ! Write the block data to a .txt file
      call open_txt_file("output/output_from_", rank, iounit)
      call print_cart_comm_type(comm_params, iounit)
      call print_block_type(ndims, u_params, iounit)
      call print_block_type(ndims, v_params, iounit)
      call print_block_type(ndims, vel_div_params, iounit)
      call print_block_type(ndims, p_params, iounit)
      call close_txt_file(iounit)

      ! Write the block data to a .dat file
      call write_block_data_to_file(u_params%data_layout, "output/u_solution.dat", comm_params%comm, u_params%matrix_ptr)
      call write_block_data_to_file(v_params%data_layout, "output/v_solution.dat", comm_params%comm, v_params%matrix_ptr)
      call write_block_data_to_file(vel_div_params%data_layout, "output/vel_div_solution.dat", &
         comm_params%comm, vel_div_params%matrix_ptr)
      call write_block_data_to_file(p_params%data_layout, "output/p_solution.dat", comm_params%comm, p_params%matrix_ptr)

      ! Deallocate used memory

      call deallocate_cart_comm_type(comm_params)

      call deallocate_finite_difference_stencil(FDstencil_params)

      call deallocate_block_type(u_params)

      call deallocate_block_type(v_params)

      call deallocate_block_type(vel_div_params)

      call deallocate_block_type(p_params)

   end subroutine TravelingWave_Poisson_2D_main

   recursive subroutine v_cycle(solver, comm, u_block, v_block, p_block, vel_div_block, FDstencil, max_level, level)
      type(SolverParamsType), intent(in) :: solver
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: u_block, v_block, p_block, vel_div_block
      type(FDstencil_type), intent(inout) :: FDstencil
      integer, intent(in) :: max_level, level

      type(block_type) :: u_block_coarse, v_block_coarse, p_block_coarse, vel_div_block_coarse
      type(FDstencil_type) :: FDstencil_coarse

      ! If max_level is reached, solve the system directly for FMG.
      if(level == max_level) then
         call solve_pressure_poisson(comm, p_block, FDstencil, vel_div_block, solver)
      else
         call solve_pressure_poisson(comm, p_block, FDstencil, vel_div_block, solver)

         call full_weighing_restriction_2D(p_block%)



      end subroutine v_cycle

      subroutine Write_analytical_function(u_block, v_block)
         type(block_type), intent(inout) :: u_block, v_block

         integer :: ii, jj
         integer, dimension(u_block%ndims) :: local_indices, global_indices
         real, dimension(u_block%ndims) :: point
         real :: x, sigma
         real :: uu, ww, eta, etax, etaxx, pp

         !$omp parallel do collapse(2) default(none) &
         !$omp shared(u_block, v_block) &
         !$omp private(ii, jj, local_indices, global_indices, point, x, sigma, uu, ww, eta, etax, etaxx, pp)
         do ii = u_block%block_begin_c(2)+1, u_block%block_end_c(2)
            do jj = u_block%block_begin_c(1)+1, u_block%block_end_c(1)
               local_indices = [jj,ii]
               global_indices = u_block%extended_global_begin_c + local_indices

               call calculate_point(u_block%ndims, -u_block%global_grid_begin, global_indices, &
                  u_block%domain_begin, u_block%grid_dx, point)

               x = point(1)
               sigma = point(2) ! we should do 1-point(2) to start from the surface! Also do it in python code

               call TravelingWave1D(Hwave, cwave, kwave, sigma, hd, wwave, start_t, x, uu, ww, eta, etax, etaxx, pp)

               u_block%matrix_ptr_2D(jj,ii) = uu
               v_block%matrix_ptr_2D(jj,ii) = ww

            end do
         end do

         !$omp end parallel do

      end subroutine Write_analytical_function

      ! Calculate velocity divergence
      subroutine calculate_velocity_divergence(u_block, v_block, FDstencil, vel_div_block)
         type(block_type), intent(in) :: u_block, v_block
         type(FDstencil_type), target, intent(inout) :: FDstencil
         type(block_type), intent(inout) :: vel_div_block

         integer :: ii, jj
         integer, dimension(u_block%ndims) :: local_indices, alpha, beta
         real, contiguous, dimension(:), pointer :: coefficients, dx, ds
         real, contiguous, dimension(:,:), pointer :: dx_2D, ds_2D
         real :: u_x, v_s

         call calculate_scaled_coefficients(u_block%ndims, u_block%extended_grid_dx, FDstencil)

         !$omp parallel do collapse(2) default(none) &
         !$omp shared(u_block, v_block, FDstencil, vel_div_block) &
         !$omp private(ii, jj, local_indices, alpha, beta, coefficients, dx, ds, dx_2D, ds_2D, u_x, v_s)
         do ii = u_block%block_begin_c(2)+1, u_block%block_end_c(2)
            do jj = u_block%block_begin_c(1)+1, u_block%block_end_c(1)
               local_indices = [jj,ii]

               call get_FD_coefficients_from_index(u_block%ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
                  [1,1], u_block%extended_block_dims, local_indices, FDstencil%scaled_stencil_coefficients, alpha, beta, &
                  coefficients)

               dx => coefficients(1:FDstencil%num_stencil_elements)
               ds => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

               call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
               call reshape_real_1D_to_2D(FDstencil%stencil_sizes, ds, ds_2D)

               call apply_stencil_2D(dx_2D, u_block%matrix_ptr_2D, local_indices, alpha, beta, u_x)
               call apply_stencil_2D(ds_2D, v_block%matrix_ptr_2D, local_indices, alpha, beta, v_s)

               new_val = u_x/Lx + v_s/hd

               vel_div_block%matrix_ptr_2D(jj,ii) = new_val
               vel_div_block%residual_matrix_ptr_2D(jj,ii) = new_val - 0.0

            end do
         end do

         !$omp end parallel do

      end subroutine calculate_velocity_divergence

      !> Write the boundary conditions for the pressure. Pressure is set to zero at the top and P_y = 0 at the bottom.
      !! Periodic boundary conditions are used for the left and right walls.
      subroutine write_pressure_BC(dims, matrix)
         integer, dimension(:), intent(in) :: dims
         real, dimension(:,:), intent(inout) :: matrix

         ! Periodic is built into the matrix, so we don't need to do anything here.


         ! At seabottom: p_y = 0 (Neumann condition)
         matrix(:,1) = matrix(:,3)


         ! At the free surface: p = 0 (Dirichlet condition) ! Replacing with a neumann condition improves the error in this part!
         !matrix(:,dims(2)) = 0.0
         matrix(:,dims(2)) = matrix(:,dims(2)-2)

      end subroutine write_pressure_BC

      subroutine solve_pressure_poisson(comm, p_block, FDstencil, vel_div_block, solver)
         type(comm_type), intent(in) :: comm
         type(block_type), intent(inout) :: p_block, vel_div_block
         type(FDstencil_type), intent(inout) :: FDstencil
         type(SolverParamsType), intent(in) :: solver

         integer :: it, converged
         real, dimension(4) :: norm_array

         call calculate_scaled_coefficients(p_block%ndims, p_block%extended_grid_dx, FDstencil)
         call write_pressure_BC(p_block%extended_block_dims, p_block%matrix_ptr_2D)

         norm_array = [1e3,1e6,1e9,1e12]

         converged = 0
         it = 0
         do while(converged /= 1 .and. it < solver%max_iter)
            call poisson_iteration(p_block, FDstencil, vel_div_block, norm_array(1))

            call sendrecv_data_neighbors(comm%comm, p_block, p_block%matrix_ptr)
            call write_pressure_BC(p_block%extended_block_dims, p_block%matrix_ptr_2D)

            call check_convergence(comm%comm, solver%tol, solver%divergence_tol, &
               1.0/product(p_block%extended_grid_size), norm_array, converged)
            if(converged == -1 .and. it > 0) then
               write(*,*) "Convergence failed"
               exit
            end if

            it = it + 1

         end do

         write(*,*) "Converged in ", it, " iterations"

      end subroutine solve_pressure_poisson

      subroutine poisson_iteration(p_block, FDstencil, vel_div_block, local_norm)
         type(block_type), intent(inout) :: p_block, vel_div_block
         type(FDstencil_type), intent(inout) :: FDstencil
         real, intent(out) :: local_norm

         integer :: ii, jj
         integer, dimension(p_block%ndims) :: p_local_indices, uv_local_indices, alpha, beta
         real, contiguous, dimension(:), pointer :: coefficients, dxx, dss
         real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
         real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
         real :: old_val, f_val, new_val

         local_norm = 0.0

         !$omp parallel do collapse(2) default(none) reduction(+:local_norm) &
         !$omp shared(p_block, vel_div_block, FDstencil) &
         !$omp private(ii, jj, p_local_indices, uv_local_indices, alpha, beta, coefficients, dxx, dss, combined_stencils, &
         !$omp combined_stencils_2D, old_val, f_val, new_val)
         do ii = p_block%block_begin_c(2)+1, p_block%block_end_c(2)
            do jj = p_block%block_begin_c(1)+1, p_block%block_end_c(1)
               p_local_indices = [jj,ii]
               uv_local_indices = vel_div_block%block_begin_c + p_local_indices - p_block%block_begin_c

               old_val = p_block%matrix_ptr_2D(jj,ii)
               f_val = vel_div_block%matrix_ptr_2D(uv_local_indices(1), uv_local_indices(2))

               call get_FD_coefficients_from_index(p_block%ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
                  [1,1], p_block%extended_block_dims, p_local_indices, &
                  FDstencil%scaled_stencil_coefficients, alpha, beta, coefficients)

               dxx => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
               dss => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

               combined_stencils = dxx/(Lx*Lx) + dss/(hd*hd)

               call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

               call update_value_from_stencil_2D(combined_stencils_2D, p_block%matrix_ptr_2D, p_local_indices, alpha, beta, &
                  f_val, new_val)

               p_block%matrix_ptr_2D(jj,ii) = new_val
               p_block%residual_matrix_ptr_2D(jj,ii) = new_val - f_val

               local_norm = local_norm + (abs(new_val - old_val)**2)

            end do
         end do

         !$omp end parallel do

      end subroutine poisson_iteration

      elemental subroutine TravelingWave1D(H, c, k, sigma, h_small, w, t, x, uu, ww, eta, etax, etaxx, pp)
         real, intent(in)  :: H, c, k, sigma, h_small, w, t, x
         real, intent(out) :: uu, ww, eta, etax, etaxx, pp

         real :: z

         z = sigma*h_small - h_small

         ! Compute analytical solution at level z
         uu    =  k*H*c/2*cosh(k*(z+h_small))/sinh(k*h_small)*cos(w*t-k*x)
         ww    = -k*H*c/2*sinh(k*(z+h_small))/sinh(k*h_small)*sin(w*t-k*x)
         eta   =  H/2*cos(w*t-k*x)   ! z = 0
         etax  =  k*H/2*sin(w*t-k*x)
         etaxx = -(k**2)*H/2*cos(w*t-k*x);
         pp    = -H*c/2*cosh(k*(z+h_small))/sinh(k*h_small)*sin(w*t-k*x); ! see p. 83 in Svendsen & Jonsson (2001)
      end subroutine TravelingWave1D

   end module TravelingWave_Poisson_2D_module
