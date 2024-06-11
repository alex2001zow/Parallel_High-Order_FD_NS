module Navier_Stokes_2D_module_old
   use constants_module, only: pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, &
      apply_FDstencil, update_value_from_stencil, calculate_scaled_coefficients, get_FD_coefficients_from_index
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type
   use functions_module, only: FunctionPtrType, FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers, &
      reshape_real_1D_to_2D
   use initialization_module, only: write_function_to_block, write_initial_condition_and_boundary

   use solver_module, only: SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, check_convergence
   implicit none

   private

   !> Number of dimensions and number of derivatives
   integer, parameter :: ndims = 2
   integer, parameter :: num_derivatives = 4
   integer, dimension(ndims*num_derivatives), parameter :: derivatives = [0,1,1,0,0,2,2,0]

   ! Number of elements in the blocks
   integer, parameter :: elements_per_index = 1, used_elements_per_index = 1

   ! Grid parameters
   integer, dimension(ndims) :: grid_size = [42,42], processor_dims = [1,1]
   logical, dimension(ndims) :: periods = [.false.,.false.]
   logical, parameter :: reorder = .true.
   real, dimension(ndims) :: domain_begin = [0,0], domain_end = [1,1]
   integer, dimension(ndims), parameter :: stencil_sizes = [3,3]
   integer, dimension(ndims), parameter :: ghost_begin = [0,0], ghost_end = [0,0]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   ! Solver parameters
   integer, parameter :: N_iterations = 500
   integer, parameter :: N_Pressure_Poisson_iterations = 50
   real, parameter :: Pressure_Poisson_tol = 1e-6, Pressure_Poisson_div_tol = 1e-1

   ! Physical parameters
   real, parameter :: system_rho = 1.0, system_nu = 0.1, system_u_lid = 1.0, system_dt = 0.001

   public :: Navier_Stokes_2D_main

contains

   !> Solve the 2D Navier-Stokes equation
   subroutine Navier_Stokes_2D_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      type(SolverParamsType):: solver_params
      type(comm_type) :: comm_params
      type(FDstencil_type) :: FDstencil_params
      type(block_type) :: u_block_params, v_block_params, p_block_params, test_block_params

      integer :: ii, iounit, converged, iter
      real, dimension(4) :: result_array_with_timings

      call set_SolverParamsType(Pressure_Poisson_tol, Pressure_Poisson_div_tol, &
         N_Pressure_Poisson_iterations, 2, solver_params)

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm_params)

      call create_block_type(ndims, elements_per_index, used_elements_per_index, &
         domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin, ghost_end, stencil_begin, stencil_end, u_block_params)

      call create_block_type(ndims, elements_per_index, used_elements_per_index, &
         domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin, ghost_end, stencil_begin, stencil_end, v_block_params)

      call create_block_type(ndims, elements_per_index, used_elements_per_index, &
         domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin, ghost_end, stencil_begin, stencil_end, p_block_params)

      !call sleeper_function(1)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      ! Scale the coefficients by dx
      call calculate_scaled_coefficients(ndims, u_block_params%extended_grid_dx, FDstencil_params)

      ! Initialize the block
      u_block_params%matrix = 0.0
      v_block_params%matrix = 0.0
      p_block_params%matrix = 0.0

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      ! Run the solver
      do ii = 1, N_iterations
         call NS_2D_one_timestep(comm_params, u_block_params, v_block_params, p_block_params, &
            FDstencil_params, solver_params, system_rho, system_nu, system_dt, converged, iter)
         write(*,*) "Iteration: ", ii, "/", N_iterations, " Converged: ", converged, " Iter: ", iter

         if(converged == -1) then
            exit
         end if

      end do

      result_array_with_timings(2) = MPI_WTIME()

      result_array_with_timings(3) = result_array_with_timings(2) - result_array_with_timings(1)

      call all_reduce_mpi_wrapper(result_array_with_timings(3), result_array_with_timings(4), 1, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_params%comm)

      ! Write out the result and timings from the master rank
      if(rank == MASTER_RANK) then
         write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(4)/world_size, " seconds"
      end if

      ! Write out our setup to a file. Just for debugging purposes
      call open_txt_file("output/output_from_", rank, iounit)

      call print_cart_comm_type(comm_params, iounit)

      call print_block_type(ndims, u_block_params, iounit)
      call print_block_type(ndims, v_block_params, iounit)
      call print_block_type(ndims, p_block_params, iounit)

      !call print_finite_difference_stencil(ndims, FDstencil_params, iounit)

      call close_txt_file(iounit)

      ! ! Write out system solution to a file
      call write_block_data_to_file(u_block_params%data_layout, "output/u_solution.dat", comm_params%comm, u_block_params%matrix)
      call write_block_data_to_file(v_block_params%data_layout, "output/v_solution.dat", comm_params%comm, v_block_params%matrix)
      call write_block_data_to_file(p_block_params%data_layout, "output/p_solution.dat", comm_params%comm, p_block_params%matrix)

      ! Deallocate data

      call deallocate_cart_comm_type(comm_params)

      call deallocate_block_type(u_block_params)

      call deallocate_block_type(v_block_params)

      call deallocate_block_type(p_block_params)

      call deallocate_finite_difference_stencil(FDstencil_params)

   end subroutine Navier_Stokes_2D_main

   ! Solve the Navier-Stokes in 2D using the fractional step method.
   subroutine NS_2D_one_timestep(comm, u_block, v_block, p_block, FDstencil, solver_params_in, rho, nu, dt, converged, it)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver_params_in
      real, intent(in) :: rho, nu, dt
      integer, intent(inout) :: converged
      integer, intent(inout) :: it

      real, dimension(product(u_block%extended_block_dims)), target :: u_rhs, v_rhs, u_star, v_star, u_x_v_y, p_temp
      real, contiguous, dimension(:,:), pointer :: u_rhs_2D, v_rhs_2D, u_star_2D, v_star_2D, u_x_v_y_2D, p_temp_2D

      real, dimension(4) :: norm_array

      call reshape_real_1D_to_2D(u_block%extended_block_dims, u_star, u_star_2D)
      call reshape_real_1D_to_2D(u_block%extended_block_dims, v_star, v_star_2D)

      norm_array = [1e12,1e18,1e36,1e48]

      u_rhs = 0
      v_rhs = 0
      u_star = 0
      v_star = 0
      u_x_v_y = 0

      ! Calculate intermediate velocity fields
      call compute_rhs(u_block%ndims, u_block%extended_local_size, u_block%extended_block_begin_c+1, used_elements_per_index, &
         u_block%matrix, v_block%matrix, FDstencil, nu, u_rhs, v_rhs)
      u_star = u_block%matrix + dt * u_rhs
      v_star = v_block%matrix + dt * v_rhs

      call write_velocity_BC(u_block%extended_block_dims, u_star_2D, v_star_2D)

      ! Calculate velocity divergence
      call calculate_velocity_divergence(u_block%ndims, u_block%extended_block_dims, u_block%block_begin_c+1, &
         used_elements_per_index, u_star, v_star, FDstencil, u_x_v_y)

      p_temp = p_block%matrix

      ! Calculate pressure correction
      call solve_poisson_problem(comm%comm, p_block%ndims, p_block%extended_block_dims, u_block%extended_block_dims, &
         p_block%extended_block_begin_c+1, used_elements_per_index, &
         p_block%matrix, p_temp, u_x_v_y, FDstencil, solver_params_in%tol, solver_params_in%max_iter, &
         norm_array, rho, dt, converged, it)

      ! Update the velocity fields and correct the pressure
      call update_velocity_fields(ndims, u_block, v_block, p_block, FDstencil, p_block%matrix, u_star, v_star, rho, dt)

      ! Write out the boundary conditions on the velocity fields
      call write_velocity_BC(u_block%extended_block_dims, u_block%matrix_ptr_2D, v_block%matrix_ptr_2D)

   end subroutine NS_2D_one_timestep

   ! Write the lid cavity boundary conditions on the velocity fields.
   subroutine write_velocity_BC(dims, u_matrix, v_matrix)
      integer, dimension(:), intent(in) :: dims
      real, dimension(:,:), intent(inout) :: u_matrix, v_matrix

      u_matrix(1,:) = 0.0
      u_matrix(:,1) = system_u_lid
      u_matrix(dims(1),:) = 0.0
      u_matrix(:,dims(2)) = 0.0

      v_matrix(1,:) = 0.0
      v_matrix(:,1) = 0.0
      v_matrix(dims(1),:) = 0.0
      v_matrix(:,dims(2)) = 0.0

   end subroutine write_velocity_BC

   ! Write the lid cavity boundary conditions on the pressure field.
   ! We want u_x = 0 at the left and right wall and u_y = 0 at the bottom wall with a Dirichlet condition on the top wall with p = 0.
   ! Not sure about the corners, perhaps u_x + u_y = 0.
   subroutine write_pressure_BC(ndims, FDstencil, p_dims, p_matrix)
      integer, intent(in) :: ndims
      type(FDstencil_type), target, intent(inout) :: FDstencil
      integer, dimension(:), intent(in) :: p_dims
      real, dimension(:), intent(inout) :: p_matrix

      integer :: ii, jj, global_index
      integer, dimension(ndims) :: p_index, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dfx, dfy
      real :: old_val, new_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(ndims, p_dims, FDstencil, p_matrix) &
      !$omp private(ii, jj, p_index, global_index, alpha, beta, coefficients, dfx, dfy, old_val, new_val)
      do ii = 1, p_dims(1)
         do jj = 1, p_dims(2)
            p_index = [ii, jj]
            call IDX_XD(ndims, p_dims, p_index, global_index)

            call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               [1, 1], p_dims, p_index, FDstencil%scaled_stencil_coefficients, alpha, beta, coefficients)

            dfx => coefficients(1:FDstencil%num_stencil_elements)
            dfy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            old_val = 0.0

            ! Left and right walls: u_x = 0 (Neumann condition)
            if (jj == 1 .or. jj == p_dims(2)) then
               call update_value_from_stencil(ndims, 1, 0, FDstencil%stencil_sizes, alpha, dfx, &
                  p_dims, p_index, p_matrix, old_val, new_val)

               p_matrix(global_index) = new_val
            end if

            ! Top and bottom walls: u_y = 0 (Neumann condition)
            if (ii == 1 .or. ii == p_dims(1)) then
               call update_value_from_stencil(ndims, 1, 0, FDstencil%stencil_sizes, alpha, dfy, &
                  p_dims, p_index, p_matrix, old_val, new_val)

               p_matrix(global_index) = new_val
            end if

            ! Top wall: p = 0 (Dirichlet condition)
            if (ii == 1) then
               p_matrix(global_index) = 0.0
            end if
         end do
      end do

      !$omp end parallel do

   end subroutine write_pressure_BC

   ! Calculate rhs for the 2D Navier-Stokes equation
   subroutine compute_rhs(ndims, dims, start_dims, num_elements, u_matrix, v_matrix, FDstencil, nu, u_rhs, v_rhs)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: dims, start_dims
      real, dimension(:), intent(in) :: u_matrix, v_matrix
      type(FDstencil_type), target, intent(inout) :: FDstencil
      real, intent(in) :: nu
      real, dimension(:), intent(inout) :: u_rhs, v_rhs

      integer :: ii, jj, global_index
      integer, dimension(ndims) :: index, alphas, betas
      real, dimension(product(FDstencil%stencil_sizes)) :: combined_stencil
      real, contiguous, dimension(:), pointer :: coefficients, dx, dy, dxx, dyy
      real, contiguous, dimension(:,:), pointer :: dx_2D, dy_2D, dxx_2D, dyy_2D
      real :: u, v, f_u, f_v

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(dims, u_matrix, v_matrix, FDstencil, nu, u_rhs, v_rhs, num_elements, ndims, start_dims) &
      !$omp private(ii, jj, index, global_index, alphas, betas, combined_stencil, coefficients, dx, dy, dxx, dyy, u, v, f_u, f_v)
      do ii = 2, dims(1)-1
         do jj = 2, dims(2)-1
            index = [ii,jj]
            call IDX_XD(ndims, dims, index, global_index)

            f_u = 0
            f_v = 0

            u = u_matrix(global_index)
            v = v_matrix(global_index)

            call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               start_dims, dims, index, FDstencil%scaled_stencil_coefficients, alphas, betas, coefficients)

            dx => coefficients(1:FDstencil%num_stencil_elements)
            dy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)
            dxx => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
            dyy => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

            !call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
            ! Find the matrix elements by using
            ! u_matrix(index(1)-alpha(1):index(1)+beta(1), index(2)-alpha(2):index(2)+beta(2)) or something similar. Get beta from get_FD_coefficients_from_index

            combined_stencil = -(u*dx + v*dy) + nu * (dxx + dyy) + f_u
            call apply_FDstencil(ndims, num_elements, 0, FDstencil%stencil_sizes, alphas, combined_stencil, dims, index, &
               u_matrix,  u_rhs(global_index))

            combined_stencil = -(u*dx + v*dy) + nu * (dxx + dyy) + f_v
            call apply_FDstencil(ndims, num_elements, 0, FDstencil%stencil_sizes, alphas, combined_stencil, dims, index, &
               v_matrix, v_rhs(global_index))

         end do
      end do

      !$omp end parallel do

   end subroutine compute_rhs

   ! Calculate velocity divergence
   subroutine calculate_velocity_divergence(ndims, dims, start_dims, num_elements, u_matrix, v_matrix, FDstencil, u_x_v_y)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: dims, start_dims
      real, dimension(:), intent(in) :: u_matrix, v_matrix
      type(FDstencil_type), target, intent(inout) :: FDstencil
      real, dimension(:), intent(inout) :: u_x_v_y

      integer :: ii, jj, global_index
      integer, dimension(ndims) :: index, alphas, betas
      integer, contiguous, dimension(:), pointer :: stencil_size
      real, contiguous, dimension(:), pointer :: coefficients, dfx, dfy
      real :: u_x, v_y

      stencil_size => FDstencil%stencil_sizes

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(dims, ndims, num_elements, start_dims, FDstencil, stencil_size, u_matrix, v_matrix, u_x_v_y) &
      !$omp private(ii, jj, index, global_index, alphas, betas, coefficients, dfx, dfy, u_x, v_y)
      do ii = 2, dims(1)-1
         do jj = 2, dims(2)-1
            index = [ii,jj]
            call IDX_XD(ndims, dims, index, global_index)

            call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               start_dims, dims, index, FDstencil%scaled_stencil_coefficients, alphas, betas, coefficients)

            dfx => coefficients(1:FDstencil%num_stencil_elements)
            dfy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfx, dims, index, u_matrix, u_x)
            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfy, dims, index, v_matrix, v_y)

            u_x_v_y(global_index) = u_x + v_y

         end do
      end do

      !$omp end parallel do

   end subroutine calculate_velocity_divergence

   subroutine solve_poisson_problem(comm, ndims, p_dims, dims, p_start_dims, num_elements, p_matrix, p_temp_matrix, u_x_v_y, &
      FDstencil, tol, max_iter, norm_array, rho, dt, converged, it)
      integer, intent(in) :: comm, ndims, num_elements, max_iter
      integer, dimension(ndims), intent(in) :: p_dims, dims, p_start_dims
      real, dimension(:), target, intent(inout) :: p_matrix, p_temp_matrix
      real, dimension(:), intent(inout) :: norm_array, u_x_v_y
      type(FDstencil_type), target, intent(inout) :: FDstencil
      real, intent(in) :: tol, rho, dt
      integer, intent(out) :: converged, it

      integer :: ii, jj, global_index
      integer, dimension(ndims) :: p_index, alphas, betas
      integer, dimension(:), pointer :: stencil_size

      real, contiguous, dimension(:), pointer :: coefficients, dfxx, dfyy, p_ptr, p_temp_ptr
      real, dimension(product(FDstencil%stencil_sizes)) :: combined_stencil
      real :: local_norm, old_val, new_val, f_val

      stencil_size => FDstencil%stencil_sizes

      p_ptr => p_matrix
      p_temp_ptr => p_temp_matrix

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < max_iter)

         local_norm = 0.0

         !$omp parallel do collapse(2) default(none) reduction(+:local_norm) &
         !$omp shared(p_dims, ndims, rho, dt, u_x_v_y, p_ptr, p_temp_ptr, FDstencil, stencil_size, num_elements, p_start_dims) &
         !$omp private(ii, jj, p_index, global_index, old_val, alphas, betas, coefficients, dfxx, dfyy, &
         !$omp combined_stencil, f_val, new_val)
         do ii = 2, p_dims(1)-1
            do jj = 2, p_dims(2)-1
               p_index = [ii,jj]
               call IDX_XD(ndims, p_dims, p_index, global_index)

               old_val = p_ptr(global_index)

               call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
                  p_start_dims, p_dims, p_index, FDstencil%scaled_stencil_coefficients, alphas, betas, coefficients)

               dfxx => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
               dfyy => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)
               combined_stencil = (dfxx + dfyy)

               f_val = (rho/dt) * u_x_v_y(global_index)

               call update_value_from_stencil(ndims, num_elements, 0, stencil_size, alphas, combined_stencil, &
                  p_dims, p_index, p_ptr, f_val, new_val)

               p_temp_ptr(global_index) = new_val

               local_norm = local_norm + (abs(new_val - old_val))**2

            end do
         end do

         !$omp end parallel do

         norm_array(1) = local_norm

         call write_pressure_BC(ndims, FDstencil, p_dims, p_ptr)

         call check_convergence(comm, tol, 1e6, 1.0/product(p_dims), norm_array, converged)
         if(converged == -1) then
            exit
         end if

         call swap_pointers(p_ptr, p_temp_ptr)

         it = it + 1

      end do

   end subroutine solve_poisson_problem

   subroutine update_velocity_fields(ndims, u_block, v_block, p_block, FDstencil, p_corr_array, u_star, v_star, rho, dt)
      integer, intent(in) :: ndims
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      real, dimension(:), intent(in) :: p_corr_array, u_star, v_star
      real, intent(in) :: rho, dt

      integer :: ii, jj, global_index, p_global_index, num_elements
      integer, dimension(ndims) :: p_start_dims, p_index, alphas, betas
      integer, contiguous, dimension(:), pointer :: stencil_size
      real, contiguous, dimension(:), pointer :: coefficients, dfx, dfy
      real :: p_corr_x, p_corr_y

      stencil_size => FDstencil%stencil_sizes
      p_start_dims = p_block%block_begin_c+1

      ! Update the velocity fields

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(p_block, ndims, num_elements, p_start_dims, FDstencil, stencil_size, &
      !$omp p_corr_array, u_block, u_star, v_block, v_star, dt, rho) &
      !$omp private(ii, jj, p_index, p_global_index, alphas, betas, coefficients, dfx, dfy, p_corr_x, p_corr_y)
      do ii = 1,p_block%extended_local_size(1)
         do jj = 1,p_block%extended_local_size(2)
            p_index = [ii,jj]
            call IDX_XD(ndims, p_block%extended_local_size, p_index, p_global_index)

            call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               p_start_dims, p_block%extended_local_size, p_index, FDstencil%scaled_stencil_coefficients, alphas, betas, &
               coefficients)

            dfx => coefficients(1:FDstencil%num_stencil_elements)
            dfy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfx, &
               p_block%extended_local_size, p_index, p_corr_array, p_corr_x)
            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfy, &
               p_block%extended_local_size, p_index, p_corr_array, p_corr_y)

            u_block%matrix(p_global_index) = u_star(p_global_index) - (dt/rho) * p_corr_x
            v_block%matrix(p_global_index) = v_star(p_global_index) - (dt/rho) * p_corr_y

         end do
      end do

      !$omp end parallel do

   end subroutine update_velocity_fields

end module Navier_Stokes_2D_module_old
