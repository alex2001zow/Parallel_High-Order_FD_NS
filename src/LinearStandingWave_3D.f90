module LinearStandingWave_3D_module
   use constants_module, only: pi
   use utility_functions_module, only: IDX_XD, IDX_XD_INV, open_txt_file, close_txt_file, sleeper_function
   use mpi_wrapper_module, only: write_block_data_to_file
   use solver_module, only: SolverParamsType, set_SolverParamsType
   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, print_cart_comm_type
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type, &
      sendrecv_data_neighbors
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      calculate_scaled_coefficients, get_FD_coefficients_from_index, update_value_from_stencil, apply_FDstencil
   use functions_module, only: calculate_point
   use solver_module, only: check_convergence

   implicit none

   private

   !> Physical parameters
   real, parameter :: g = 9.81, mu = 1.3059*(1e-6), rho = 1000.0, nu = mu/rho

   !> Domain parameters
   real, parameter :: Ls = 1.0, Ly = 31.0, Lx = 31.0

   !> Wave parameters
   real, parameter :: start_t = 1.0
   real, parameter :: k_y = 2*pi/Ly, k_x = 2*pi/Lx, k = sqrt(k_y**2 + k_x**2), hs = 1, c = sqrt(g/k*tanh(k*hs)), &
      omega = c*k, H = 2

   !> Grid parameters
   integer, parameter :: ndims = 3
   integer, dimension(ndims), parameter :: grid_size = [8,8,8], processor_dims = [1,1,1]
   logical, dimension(ndims), parameter :: periods = [.false., .false., .true.]
   logical, parameter :: reorder = .true.
   real, dimension(ndims), parameter :: domain_begin = [0,0,0], domain_end = [Ls,Ly,Lx]
   integer, dimension(ndims), parameter :: stencil_sizes = [3,3,3]
   integer, dimension(ndims), parameter :: ghost_begin = [1,0,0], ghost_end = [1,0,0]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   !> Solver parameters
   real, parameter :: tol = 1e-10, div_tol = 1e-1
   integer, parameter :: max_iter = 1000

   public :: LinearStandingWave_3D_main

contains

   subroutine LinearStandingWave_3D_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      integer, parameter :: num_derivatives = 6
      integer, dimension(ndims*num_derivatives), parameter :: derivatives = &
         [0,0,1,& ! First derivative in x
         0,1,0,& ! First derivative in y
         1,0,0,& ! First derivative in z (sigma)
         0,0,2,& ! Second derivative in x
         0,2,0,& ! Second derivative in y
         2,0,0] ! Second derivative in z (sigma)

      type(SolverParamsType) :: solver_params
      type(comm_type) :: comm_params
      type(block_type) :: w_params, v_params, u_params, vel_div_params, p_params
      type(FDstencil_type) :: FDstencil_params

      integer :: iounit

      !call sleeper_function(1)

      ! Set the solver parameters: tol, div_tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(tol, div_tol, max_iter, 2, solver_params)

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin*0, ghost_end*0, stencil_begin, stencil_end, w_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin*0, ghost_end*0, stencil_begin, stencil_end, v_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin*0, ghost_end*0, stencil_begin, stencil_end, u_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin*0, ghost_end*0, stencil_begin, stencil_end, vel_div_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin, ghost_end, stencil_begin, stencil_end, p_params)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      w_params%matrix_ptr = 0.0
      v_params%matrix_ptr = 0.0
      u_params%matrix_ptr = 0.0
      vel_div_params%matrix_ptr = 0.0
      p_params%matrix_ptr = 0.0

      call Write_analytical_function(w_params, v_params, u_params)

      call sendrecv_data_neighbors(comm_params%comm, w_params, w_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, v_params, v_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, u_params, u_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, p_params, p_params%matrix_ptr)

      call calculate_velocity_divergence(w_params, v_params, u_params, FDstencil_params, vel_div_params)

      ! Time and solve the system

      call solve_pressure_poisson(comm_params, p_params, FDstencil_params, vel_div_params, solver_params)

      call sendrecv_data_neighbors(comm_params%comm, w_params, w_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, v_params, v_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, u_params, u_params%matrix_ptr)
      call sendrecv_data_neighbors(comm_params%comm, p_params, p_params%matrix_ptr)

      ! Write the block data to a .txt file
      call open_txt_file("output/output_from_", rank, iounit)
      call print_cart_comm_type(comm_params, iounit)
      call print_block_type(ndims, w_params, iounit)
      call print_block_type(ndims, v_params, iounit)
      call print_block_type(ndims, u_params, iounit)
      call print_block_type(ndims, vel_div_params, iounit)
      call print_block_type(ndims, p_params, iounit)
      call close_txt_file(iounit)

      ! Write the block data to a .dat file
      call write_block_data_to_file(w_params%data_layout, "output/w_solution.dat", comm_params%comm, w_params%matrix_ptr)
      call write_block_data_to_file(v_params%data_layout, "output/v_solution.dat", comm_params%comm, v_params%matrix_ptr)
      call write_block_data_to_file(u_params%data_layout, "output/u_solution.dat", comm_params%comm, u_params%matrix_ptr)
      call write_block_data_to_file(vel_div_params%data_layout, "output/vel_div_solution.dat", &
         comm_params%comm, vel_div_params%matrix_ptr)
      call write_block_data_to_file(p_params%data_layout, "output/p_solution.dat", comm_params%comm, p_params%matrix_ptr)

      ! Deallocate used memory

      call deallocate_block_type(w_params)
      call deallocate_block_type(v_params)
      call deallocate_block_type(u_params)
      call deallocate_block_type(vel_div_params)
      call deallocate_block_type(p_params)

      call deallocate_cart_comm_type(comm_params)

      call deallocate_finite_difference_stencil(FDstencil_params)

   end subroutine LinearStandingWave_3D_main

   subroutine solve_pressure_poisson(comm, p_block, FDstencil, vel_div_block, solver)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: p_block, vel_div_block
      type(FDstencil_type), intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver

      integer :: it, converged
      real, dimension(4) :: norm_array

      call calculate_scaled_coefficients(p_block%ndims, p_block%extended_grid_dx, FDstencil)
      call write_pressure_BC(p_block%ndims, p_block%extended_block_dims, p_block%matrix_ptr)

      norm_array = [1e3,1e6,1e9,1e12]

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < solver%max_iter)
         call poisson_iteration(p_block, FDstencil, vel_div_block, norm_array(1))

         call sendrecv_data_neighbors(comm%comm, p_block, p_block%matrix_ptr)
         call write_pressure_BC(p_block%ndims, p_block%extended_block_dims, p_block%matrix_ptr)

         call check_convergence(comm%comm, solver%tol, solver%divergence_tol, &
            1.0/product(p_block%extended_grid_size), norm_array, converged)
         if(converged == -1 .and. it > 0) then
            write(*,*) "Convergence failed"
            exit
         end if

         it = it + 1
         converged = 0

      end do

      write(*,*) "Converged in ", it, " iterations"

   end subroutine solve_pressure_poisson

   subroutine poisson_iteration(p_block, FDstencil, vel_div_block, local_norm)
      type(block_type), intent(inout) :: p_block, vel_div_block
      type(FDstencil_type), intent(inout) :: FDstencil
      real, intent(out) :: local_norm

      integer :: ii, jj, kk, p_global_index, uv_global_index
      integer, dimension(p_block%ndims) :: p_local_indices, uv_local_indices, alphas
      real, dimension(:), pointer :: coefficients, dxx, dyy, dss
      real, dimension(product(FDstencil%stencil_sizes)) :: combined_stencils
      real :: old_val, f_val, new_val

      local_norm = 0.0
      do ii = p_block%block_begin_c(1)+1, p_block%block_end_c(1)
         do jj = p_block%block_begin_c(2)+1, p_block%block_end_c(2)
            do kk = p_block%block_begin_c(3)+1, p_block%block_end_c(3)
               p_local_indices = [ii,jj]
               uv_local_indices = vel_div_block%block_begin_c + p_local_indices - p_block%block_begin_c
               call IDX_XD(p_block%ndims, p_block%extended_block_dims, p_local_indices, p_global_index)
               call IDX_XD(vel_div_block%ndims, vel_div_block%extended_block_dims, uv_local_indices, uv_global_index)
               !print *, p_local_indices, uv_local_indices

               old_val = p_block%matrix_ptr(p_global_index)

               call get_FD_coefficients_from_index(p_block%ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
                  p_block%extended_block_begin_c+1, p_block%extended_block_dims, p_local_indices, &
                  FDstencil%scaled_stencil_coefficients, alphas, coefficients)

               dxx => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)
               dyy => coefficients(4 * FDstencil%num_stencil_elements + 1:5 * FDstencil%num_stencil_elements)
               dss => coefficients(5 * FDstencil%num_stencil_elements + 1:6 * FDstencil%num_stencil_elements)

               combined_stencils = dxx + dyy + (1.0/(H**2))*dss

               f_val = vel_div_block%matrix_ptr(uv_global_index)

               call update_value_from_stencil(p_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, combined_stencils, &
                  p_block%extended_block_dims, p_local_indices, p_block%matrix_ptr, f_val, new_val)

               p_block%matrix_ptr(p_global_index) = new_val

               local_norm = local_norm + (abs(new_val - old_val)**2)

            end do
         end do
      end do

   end subroutine poisson_iteration

   ! Calculate velocity divergence
   subroutine calculate_velocity_divergence(w_block, v_block, u_block, FDstencil, vel_div_block)
      type(block_type), intent(in) :: w_block, v_block, u_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(block_type), intent(inout) :: vel_div_block

      integer :: ii, jj, kk, vel_global_index
      integer, dimension(u_block%ndims) :: indices, alphas
      real, dimension(:), pointer :: coefficients, dx, dy, ds
      real :: u_x, v_y, w_s

      call calculate_scaled_coefficients(u_block%ndims, u_block%extended_grid_dx, FDstencil)

      !$omp parallel do collapse(3) default(none) &
      !$omp shared(w_block, v_block, u_block, FDstencil, vel_div_block) &
      !$omp private(ii, jj, kk, indices, vel_global_index, alphas, coefficients, dx, dy, ds, u_x, v_y, w_s)
      do ii = u_block%block_begin_c(1)+1, u_block%block_end_c(1)
         do jj = u_block%block_begin_c(2)+1, u_block%block_end_c(2)
            do kk = u_block%block_begin_c(3)+1, u_block%block_end_c(3)
               indices = [ii,jj,kk]
               call IDX_XD(u_block%ndims, u_block%extended_block_dims, indices, vel_global_index)

               call get_FD_coefficients_from_index(u_block%ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
                  [1,1,1], u_block%extended_block_dims, indices, FDstencil%scaled_stencil_coefficients, alphas, coefficients)

               dx => coefficients(1:FDstencil%num_stencil_elements)
               dy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)
               ds => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)

               call apply_FDstencil(u_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, dx, &
                  u_block%extended_block_dims, indices, u_block%matrix_ptr, u_x)
               call apply_FDstencil(u_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, dy, &
                  u_block%extended_block_dims, indices, v_block%matrix_ptr, v_y)
               call apply_FDstencil(u_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, ds, &
                  u_block%extended_block_dims, indices, w_block%matrix_ptr, w_s)

               vel_div_block%matrix_ptr(vel_global_index) = u_x + v_y + (1.0/H)*w_s

            end do
         end do
      end do

      !$omp end parallel do

   end subroutine calculate_velocity_divergence

   !> Write the boundary conditions for the pressure. Pressure is set to zero at the top and P_y = 0 at the bottom.
   !! Periodic boundary conditions are used for the left and right walls.
   subroutine write_pressure_BC(ndims, dims, matrix)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: dims
      real, dimension(:), intent(inout) :: matrix

      integer :: ii, jj, kk, global_index, new_global_index
      integer, dimension(ndims) :: indices

      !$omp parallel do collapse(3) default(none) &
      !$omp shared(ndims, dims, matrix) &
      !$omp private(ii, jj, kk, indices, global_index, new_global_index)
      do ii = 1, dims(1)
         do jj = 1, dims(2)
            do kk = 1, dims(3)
               indices = [ii, jj, kk]
               call IDX_XD(ndims, dims, indices, global_index)

               ! Top wall: p = 0 (Dirichlet condition)
               if (ii == 1) then
                  matrix(global_index) = 0.0
               end if

               ! Bottom wall: p_y = 0 (Neumann condition)
               if (ii == dims(1)) then
                  call IDX_XD(ndims, dims, [ii-2, jj, kk], new_global_index)

                  matrix(global_index) = matrix(new_global_index)
               end if

            end do
         end do
      end do

      !$omp end parallel do

   end subroutine write_pressure_BC

   subroutine Write_analytical_function(w_block, v_block, u_block)
      type(block_type), intent(inout) :: w_block, v_block, u_block

      integer :: ii, jj, kk, global_index
      integer, dimension(u_block%ndims) :: local_indices, global_indices
      real, dimension(u_block%ndims) :: point
      real :: sigma, y, x
      real :: eta, w, v, u

      do ii = u_block%block_begin_c(1)+1, u_block%block_end_c(1)
         do jj = u_block%block_begin_c(2)+1, u_block%block_end_c(2)
            do kk = u_block%block_begin_c(3)+1, u_block%block_end_c(3)
               local_indices = [ii,jj,kk]
               global_indices = u_block%global_begin_c + local_indices
               call IDX_XD(u_block%ndims, u_block%extended_block_dims, local_indices, global_index)

               call calculate_point(u_block%ndims, [0,0,0], global_indices, &
                  u_block%domain_begin, u_block%grid_dx, point)
               !print *, global_index, local_indices, global_indices, point

               sigma = point(1)
               y = point(2)
               x = point(3)

               call LinearStandingWave_2D(sigma, y, x, start_t, k_y, k_x, k, hs, c, omega, H, eta, w, v, u)

               w_block%matrix_ptr(global_index) = w
               v_block%matrix_ptr(global_index) = v
               u_block%matrix_ptr(global_index) = u

            end do
         end do
      end do

   end subroutine Write_analytical_function

   !> Calculate the analytical solution for the 2D travelling wave. This being the free surface and the velocity field.
   elemental subroutine LinearStandingWave_2D(sigma, y, x, t, k_y, k_x, k, hs, c, omega, H, eta, w, v, u)
      real, intent(in) :: sigma, y, x, t
      real, intent(in) :: k_y, k_x, k, hs, c, omega, H
      real, intent(out) :: eta, w, v, u

      real :: z

      ! Compute the level z from sigma [0,1]
      z = (sigma - 1.0)* hs

      ! Compute analytical free surface. Always at z = 0
      eta = H / 2 * cos(omega * t) * cos(k_y * y) * cos(k_x * x)

      ! Compute the analytical velocity solution at level z. Found from phi(x,y,z,t) = -(H/2) * c * cosh(k*(z+hs)) * sin(omega*t) * cos(k_y*y) * cos(k_x*x)/sinh(k*hs)
      w = -H * c * k * sin(omega * t) * cos(k_y * y) * cos(k_x * x) * sinh(k * (z + hs)) / (2 * sinh(k * hs))
      v = H * c * k_y * sin(k_y * y) * sin(omega * t) * cos(k_x * x) * cosh(k * (z + hs)) / (2 * sinh(k * hs))
      u = H * c * k_x * sin(k_x * x) * sin(omega * t) * cos(k_y * y) * cosh(k * (z + hs)) / (2 * sinh(k * hs))

   end subroutine LinearStandingWave_2D

end module LinearStandingWave_3D_module

