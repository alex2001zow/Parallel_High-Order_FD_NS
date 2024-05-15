module Navier_Stokes_2D_module
   use constants_module, only: pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type, print_cartesian_grid
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, &
      apply_FDstencil, update_value_from_stencil, calculate_scaled_coefficients, get_FD_coefficients_from_index
   use block_module, only: block_type, create_block_type, deallocate_block_type, sendrecv_data_neighbors, print_block_type
   use functions_module, only: FunctionPtrType, FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers
   use initialization_module, only: write_function_to_block, write_initial_condition_and_boundary

   use solver_module, only: SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, check_convergence
   implicit none

   private

   integer, parameter :: N_iterations = 20
   integer, parameter :: N_Pressure_Poisson_iterations = 50
   real, parameter :: Pressure_Poisson_tol = 1e-6, Pressure_Poisson_divergence_tol = 1e-1
   real, parameter :: system_rho = 1.0, system_nu = 0.1, system_u_lid = 1.0, system_dt = 0.0001
   public :: Navier_Stokes_2D_main

contains

   subroutine Navier_Stokes_2D_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      integer :: ndims
      integer, dimension(:), allocatable :: grid_size, processor_dims, stencil_sizes
      real, dimension(:), allocatable :: domain_begin, domain_end
      type(SolverParamsType) :: solver_params

      ndims = 2

      allocate(grid_size(ndims))
      allocate(processor_dims(ndims))
      allocate(domain_begin(ndims))
      allocate(domain_end(ndims))
      allocate(stencil_sizes(ndims))

      grid_size = 42
      processor_dims = 1
      domain_begin = 0
      domain_end = 1
      stencil_sizes = 3

      ! Set the solver parameters tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(Pressure_Poisson_tol, Pressure_Poisson_divergence_tol, N_Pressure_Poisson_iterations, &
         2, solver_params)

      call Navier_Stokes_2D_analytical(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         stencil_sizes, solver_params)

   end subroutine Navier_Stokes_2D_main

   !> The analytical solution of the 2D Navier-Stokes equation
   subroutine Navier_Stokes_2D_analytical(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      stencil_sizes, solver_params)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(2), intent(in) :: grid_size, processor_dims, stencil_sizes
      real, dimension(2), intent(in) :: domain_begin, domain_end
      type(SolverParamsType), intent(in) :: solver_params

      integer, parameter :: num_derivatives = 4
      integer, dimension(2*num_derivatives), parameter :: derivatives = [0,1,1,0,0,2,2,0]

      call Navier_Stokes_2D_analytical_wrapper(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, stencil_sizes, solver_params)

   end subroutine Navier_Stokes_2D_analytical

   !> Solve the 2D Navier-Stokes equation
   subroutine Navier_Stokes_2D_analytical_wrapper(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      num_derivatives, derivatives, stencil_sizes, solver_params)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(:), intent(in) :: grid_size, processor_dims
      real, dimension(:), intent(in) :: domain_begin, domain_end

      integer, intent(in) :: num_derivatives
      integer, dimension(:), intent(in) :: derivatives, stencil_sizes
      type(SolverParamsType), intent(in) :: solver_params

      integer :: iounit, ii, converged, iter
      real, dimension(4) :: result_array_with_timings

      real :: dt

      type(comm_type) :: comm_params
      type(FDstencil_type) :: uv_FDstencil, p_FDstencil
      type(block_type) :: u_block_params, v_block_params, p_block_params

      integer, dimension(ndims) ::  bc_begin, bc_end, ghost_begin, ghost_end, stencil_begin, stencil_end

      bc_begin = 0
      bc_end = 0
      ghost_begin = stencil_sizes/2
      ghost_end = stencil_sizes/2
      stencil_begin = stencil_sizes/2
      stencil_end = stencil_sizes/2

      call create_cart_comm_type(ndims, processor_dims, rank, world_size, comm_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         bc_begin, bc_end, ghost_begin*0, ghost_end*0, stencil_begin, stencil_end, u_block_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         bc_begin, bc_end, ghost_begin*0, ghost_end*0, stencil_begin, stencil_end, v_block_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         bc_begin*0, bc_end*0, ghost_begin, ghost_end, stencil_begin, stencil_end, p_block_params)

      !call sleeper_function(1)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, uv_FDstencil)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, p_FDstencil)

      ! Scale the coefficients by dx
      call calculate_scaled_coefficients(ndims, u_block_params%dx, uv_FDstencil)

      call calculate_scaled_coefficients(ndims, p_block_params%dx, p_FDstencil)

      allocate(p_block_params%temp_matrix_ptr(p_block_params%extended_num_elements))

      ! Initialize the block
      u_block_params%matrix_ptr = 0.0
      v_block_params%matrix_ptr = 0.0
      p_block_params%matrix_ptr = 0.0
      p_block_params%temp_matrix_ptr = 0.0

      !call test_pressure_bc(ndims, p_block_params%extended_block_dims, p_block_params%extended_block_begin_c, &
      !p_block_params%matrix_ptr)

      !call write_pressure_BC(ndims, p_block_params%extended_block_dims, p_block_params%extended_block_begin_c, &
      !p_block_params%matrix_ptr)

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      dt = system_dt
      do ii = 1, N_iterations
         call NS_2D_one_timestep(comm_params, u_block_params, v_block_params, p_block_params, &
            uv_FDstencil, p_FDstencil, solver_params, system_rho, system_nu, dt, converged, iter)
         write(*,*) "Iteration: ", ii, "/", N_iterations, " Converged: ", converged, " Iter: ", iter

         if(converged == -1) then
            !exit
         end if

      end do

      result_array_with_timings(2) = MPI_WTIME()

      result_array_with_timings(3) = result_array_with_timings(2) - result_array_with_timings(1)

      call all_reduce_mpi_wrapper(result_array_with_timings(3), result_array_with_timings(4), 1, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_params%comm)

      ! Write out the result and timings from the master rank
      if(rank == MASTER_RANK) then
         !call print_resultType(result)
         write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(4)/world_size, " seconds"
      end if

      ! Write out our setup to a file. Just for debugging purposes
      call open_txt_file("output/output_from_", rank, iounit)

      !call print_cart_comm_type(comm_params, iounit)

      !call print_finite_difference_stencil(ndims, FDstencil_params, iounit)

      call print_block_type(ndims, u_block_params, iounit)
      call print_block_type(ndims, v_block_params, iounit)
      call print_block_type(ndims, p_block_params, iounit)

      call close_txt_file(iounit)

      ! ! Write out system solution to a file
      call write_block_data_to_file(u_block_params%data_layout, "output/u_solution.dat", &
         comm_params%comm, u_block_params%matrix_ptr)
      call write_block_data_to_file(v_block_params%data_layout, "output/v_solution.dat", &
         comm_params%comm, v_block_params%matrix_ptr)
      call write_block_data_to_file(p_block_params%data_layout, "output/p_solution.dat", &
         comm_params%comm, p_block_params%matrix_ptr)

      ! Deallocate data

      call deallocate_cart_comm_type(comm_params)

      call deallocate_block_type(u_block_params)

      call deallocate_block_type(v_block_params)

      call deallocate_block_type(p_block_params)

      call deallocate_finite_difference_stencil(uv_FDstencil)

      call deallocate_finite_difference_stencil(p_FDstencil)

   end subroutine Navier_Stokes_2D_analytical_wrapper

   ! Solve the Navier-Stokes in 2D using the fractional step method.
   subroutine NS_2D_one_timestep(comm, u_block, v_block, p_block, uv_FDstencil, p_FDstencil, solver_params_in, &
      rho, nu, dt, converged, it)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: uv_FDstencil, p_FDstencil
      type(SolverParamsType), intent(in) :: solver_params_in
      real, intent(in) :: rho, nu, dt
      integer, intent(out) :: converged
      integer, intent(out) :: it

      real, dimension(product(u_block%extended_local_size)) :: u_star, v_star, u_x_v_y
      real :: u_error, v_error
      real, dimension(2) :: local_error_array, global_error_array

      integer :: ndims, num_elements

      real, dimension(4) :: norm_array

      ndims = comm%ndims
      num_elements = 1
      norm_array = [1e12,1e18,1e36,1e48]

      u_star = 0
      v_star = 0
      u_x_v_y = 0

      ! Calculate intermediate velocity fields
      call DOPRI5(u_block, v_block, uv_FDstencil, nu, dt, u_star, v_star, u_error, v_error)
      local_error_array = [u_error, v_error]

      ! Do an allreduce to get the maximum error
      call all_reduce_mpi_wrapper(local_error_array, global_error_array, 2, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm%comm)

      call sendrecv_data_neighbors(comm%comm, u_block, u_star)
      call sendrecv_data_neighbors(comm%comm, v_block, v_star)

      ! Write out the boundary conditions on the velocity fields
      call write_velocity_BC(u_block%ndims, u_block%extended_block_dims, u_block%global_begin_c, u_star, v_star)

      ! Calculate velocity divergence
      call calculate_velocity_divergence(u_block, v_block, u_star, v_star, uv_FDstencil, u_x_v_y)

      ! Calculate pressure correction. We need to have a different stencil dx here because of the ghost points.
      call solve_poisson_problem(u_block, v_block, p_block, p_block%matrix_ptr, p_block%temp_matrix_ptr, u_x_v_y, p_FDstencil, &
         solver_params_in, comm, norm_array, rho, dt, converged, it)

      ! Update the velocity fields and correct the pressure
      call update_velocity_fields(u_block, v_block, p_block, uv_FDstencil, u_star, v_star, rho, dt)

      ! Write out the boundary conditions on the velocity fields
      call write_velocity_BC(u_block%ndims, u_block%extended_block_dims, u_block%global_begin_c, &
         u_block%matrix_ptr, v_block%matrix_ptr)

      call sendrecv_data_neighbors(comm%comm, u_block, u_block%matrix_ptr)
      call sendrecv_data_neighbors(comm%comm, v_block, v_block%matrix_ptr)
      call sendrecv_data_neighbors(comm%comm, p_block, p_block%matrix_ptr)

   end subroutine NS_2D_one_timestep

   ! Write the lid cavity boundary conditions on the velocity fields.
   subroutine write_velocity_BC(ndims, dims, begin_indices, u_matrix, v_matrix)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: dims, begin_indices
      real, dimension(:), intent(inout) :: u_matrix, v_matrix

      integer :: global_index
      integer, dimension(ndims) :: indices

      !$omp parallel do default(none) &
      !$omp shared(ndims, dims, begin_indices, u_matrix, v_matrix) &
      !$omp private(global_index, indices)
      do global_index = 1, product(dims)
         call IDX_XD_INV(ndims, dims, global_index, indices)
         indices = begin_indices + indices

         if(indices(1) == 1 .or. indices(1) == dims(1) .or. indices(2) == 1 .or. indices(2) == dims(2)) then
            u_matrix(global_index) = 0.0
            v_matrix(global_index) = 0.0
         end if

         if(indices(1) == 1) then
            u_matrix(global_index) = system_u_lid
         end if
      end do

      !$omp end parallel do

   end subroutine write_velocity_BC

   ! Write the lid cavity boundary conditions on the pressure field.
   ! We want u_x = 0 at the left and right wall and u_y = 0 at the bottom wall and top wall. We anchor the pressure at the bottom corner.
   ! Not sure about the corners, perhaps u_x + u_y = 0.
   subroutine write_pressure_BC(ndims, dims, begin_indices, p_matrix)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: dims, begin_indices
      real, dimension(:), intent(inout) :: p_matrix

      integer :: global_index, boundary_index, ghost_index, interior_index
      integer, dimension(ndims) :: indices

      !$omp parallel do default(none) &
      !$omp shared(ndims, dims, begin_indices, p_matrix) &
      !$omp private(global_index, indices, boundary_index, ghost_index, interior_index)
      do global_index = 1, product(dims)
         call IDX_XD_INV(ndims, dims, global_index, indices)
         indices = begin_indices + indices ! Go from the local space to the global space

         ! Top wall: u_y = 0 (Neumann condition)
         if (indices(1) == 2) then
            boundary_index = global_index
            call IDX_XD(ndims, dims, [indices(1)-1, indices(2)], ghost_index)
            call IDX_XD(ndims, dims, [indices(1)+1, indices(2)], interior_index)

            p_matrix(ghost_index) = p_matrix(interior_index)
            p_matrix(boundary_index) = p_matrix(interior_index)
         end if

         ! Bottom wall: u_y = 0 (Neumann condition)
         if(indices(1) == dims(1)-1) then
            boundary_index = global_index
            call IDX_XD(ndims, dims, [indices(1)+1, indices(2)], ghost_index)
            call IDX_XD(ndims, dims, [indices(1)-1, indices(2)], interior_index)

            p_matrix(ghost_index) = p_matrix(interior_index)
            p_matrix(boundary_index) = p_matrix(interior_index)
         end if

         ! Left wall: u_x = 0 (Neumann condition)
         if (indices(2) == 2) then
            boundary_index = global_index
            call IDX_XD(ndims, dims, [indices(1), indices(2)-1], ghost_index)
            call IDX_XD(ndims, dims, [indices(1), indices(2)+1], interior_index)

            p_matrix(ghost_index) = p_matrix(interior_index)
            p_matrix(boundary_index) = p_matrix(interior_index)
         end if

         ! Right wall: u_x = 0 (Neumann condition)
         if(indices(2) == dims(2)-1) then
            boundary_index = global_index
            call IDX_XD(ndims, dims, [indices(1), indices(2)+1], ghost_index)
            call IDX_XD(ndims, dims, [indices(1), indices(2)-1], interior_index)

            p_matrix(ghost_index) = p_matrix(interior_index)
            p_matrix(boundary_index) = p_matrix(interior_index)
         end if

         ! Bottom corner of the true domain. Not at the corner ghost point. Unsure about this: p = 0 (Dirichlet condition)
         if (indices(1) == dims(1)-1 .and. indices(2) == dims(2)-1) then
            p_matrix(global_index) = 0.0
         end if
      end do

      !$omp end parallel do

   end subroutine write_pressure_BC

   !> Subroutine to implement the DOPRI5 method for the Navier-Stokes equation
   subroutine DOPRI5(u_block, v_block, FDstencil, nu, dt, u_star, v_star, u_error, v_error)
      type(block_type), intent(in) :: u_block, v_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      real, intent(in) :: nu, dt
      real, dimension(:), intent(inout) :: u_star, v_star
      real, intent(out) :: u_error, v_error

      real, dimension(product(u_block%extended_local_size)) :: k1_u, k1_v, k2_u, k2_v, k3_u, k3_v, k4_u, k4_v, k5_u, k5_v
      real, dimension(product(u_block%extended_local_size)) :: k6_u, k6_v, k7_u, k7_v

      ! Butcher tableau constants for Dormand-Prince method
      real, parameter :: a21 = 1.0/5.0, a31 = 3.0/40.0, a32 = 9.0/40.0
      real, parameter :: a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0
      real, parameter :: a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0, a54 = -212.0/729.0
      real, parameter :: a61 = 9017.0/3168.0, a62 = -355.0/33.0, a63 = 46732.0/5247.0, &
         a64 = 49.0/176.0, a65 = -5103.0/18656.0
      real, parameter :: a71 = 35.0/384.0, a73 = 500.0/1113.0, a74 = 125.0/192.0, &
         a75 = -2187.0/6784.0, a76 = 11.0/84.0
      real, parameter :: b1 = 5179.0/57600.0, b3 = 7571.0/16695.0, b4 = 393.0/640.0, &
         b5 = -92097.0/339200.0, b6 = 187.0/2100.0, b7 = 1.0/40.0
      real, parameter :: e1 = 35.0/384.0 - 5179.0/57600.0
      real, parameter :: e3 = 500.0/1113.0 - 7571.0/16695.0
      real, parameter :: e4 = 125.0/192.0 - 393.0/640.0
      real, parameter :: e5 = -2187.0/6784.0 + 92097.0/339200.0
      real, parameter :: e6 = 11.0/84.0 - 187.0/2100.0
      real, parameter :: e7 = -1.0/40.0

      ! Initial computation for k1 using the initial conditions or previous step values
      call compute_rhs(u_block, v_block, u_block%matrix_ptr, v_block%matrix_ptr, FDstencil, nu, k1_u, k1_v)

      ! Compute k2
      u_star = u_block%matrix_ptr + dt*a21*k1_u
      v_star = v_block%matrix_ptr + dt*a21*k1_v
      call compute_rhs(u_block, v_block, u_star, v_star, FDstencil, nu, k2_u, k2_v)

      ! Compute k3
      u_star = u_block%matrix_ptr + dt*(a31*k1_u + a32*k2_u)
      v_star = v_block%matrix_ptr + dt*(a31*k1_v + a32*k2_v)
      call compute_rhs(u_block, v_block, u_star, v_star, FDstencil, nu, k3_u, k3_v)

      ! Compute k4
      u_star = u_block%matrix_ptr + dt*(a41*k1_u + a42*k2_u + a43*k3_u)
      v_star = v_block%matrix_ptr + dt*(a41*k1_v + a42*k2_v + a43*k3_v)
      call compute_rhs(u_block, v_block, u_star, v_star, FDstencil, nu, k4_u, k4_v)

      ! Compute k5
      u_star = u_block%matrix_ptr + dt*(a51*k1_u + a52*k2_u + a53*k3_u + a54*k4_u)
      v_star = v_block%matrix_ptr + dt*(a51*k1_v + a52*k2_v + a53*k3_v + a54*k4_v)
      call compute_rhs(u_block, v_block, u_star, v_star, FDstencil, nu, k5_u, k5_v)

      ! Compute k6
      u_star = u_block%matrix_ptr + dt*(a61*k1_u + a62*k2_u + a63*k3_u + a64*k4_u + a65*k5_u)
      v_star = v_block%matrix_ptr + dt*(a61*k1_v + a62*k2_v + a63*k3_v + a64*k4_v + a65*k5_v)
      call compute_rhs(u_block, v_block, u_star, v_star, FDstencil, nu, k6_u, k6_v)

      ! Compute k7
      u_star = u_block%matrix_ptr + dt*(a71*k1_u + a73*k3_u + a74*k4_u + a75*k5_u + a76*k6_u)
      v_star = v_block%matrix_ptr + dt*(a71*k1_v + a73*k3_v + a74*k4_v + a75*k5_v + a76*k6_v)
      call compute_rhs(u_block, v_block, u_star, v_star, FDstencil, nu, k7_u, k7_v)

      ! Calculate u_star and v_star and compute the error
      u_star = u_block%matrix_ptr + dt*(b1*k1_u + b3*k3_u + b4*k4_u + b5*k5_u + b6*k6_u + b7*k7_u)
      v_star = v_block%matrix_ptr + dt*(b1*k1_v + b3*k3_v + b4*k4_v + b5*k5_v + b6*k6_v + b7*k7_v)
      u_error = sum(abs(e1*k1_u + e3*k3_u + e4*k4_u + e5*k5_u + e6*k6_u + e7*k7_u)**2)  !! This will also take the ghost points into account which is not correct.
      v_error = sum(abs(e1*k1_v + e3*k3_v + e4*k4_v + e5*k5_v + e6*k6_v + e7*k7_v)**2) !! This will also take the ghost points into account which is not correct.

   end subroutine DOPRI5

   ! Calculate rhs for the 2D Navier-Stokes equation
   subroutine compute_rhs(u_block, v_block, u_matrix, v_matrix, FDstencil, nu, u_rhs, v_rhs)
      type(block_type), intent(in) :: u_block, v_block
      real, dimension(:), intent(in) :: u_matrix, v_matrix
      type(FDstencil_type), target, intent(inout) :: FDstencil
      real, intent(in) :: nu
      real, dimension(:), intent(inout) :: u_rhs, v_rhs

      integer :: ii, jj, global_index
      integer, dimension(u_block%ndims) :: indices, alphas
      real, dimension(:), pointer :: coefficients, dfx, dfy, dfxx, dfyy
      real, dimension(product(FDstencil%stencil_sizes)) :: combined_stencils
      real :: u, v, f_u, f_v

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(u_block, v_block, u_matrix, v_matrix, FDstencil, nu, u_rhs, v_rhs) &
      !$omp private(ii, jj, indices, global_index, alphas, coefficients, dfx, dfy, dfxx, dfyy, combined_stencils, &
      !$omp u, v, f_u, f_v)
      do ii = u_block%block_begin_c(1)+2, u_block%block_end_c(1)-1
         do jj = u_block%block_begin_c(2)+2, u_block%block_end_c(2)-1
            indices = [ii,jj]
            call IDX_XD(u_block%ndims, u_block%extended_block_dims, indices, global_index)

            u = u_matrix(global_index)
            v = v_matrix(global_index)

            f_u = 0
            f_v = 0

            call get_FD_coefficients_from_index(u_block%ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               u_block%extended_block_begin_c+1, u_block%extended_block_dims, &
               indices, FDstencil%scaled_stencil_coefficients, alphas, coefficients)

            dfx => coefficients(1:FDstencil%num_stencil_elements)
            dfy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)
            dfxx => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
            dfyy => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

            combined_stencils = -(u * dfx + v * dfy) + nu * (dfxx + dfyy) + f_u

            call apply_FDstencil(u_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, &
               combined_stencils, u_block%extended_block_dims, indices, u_matrix, u_rhs(global_index))

            combined_stencils = -(u * dfx + v * dfy) + nu * (dfxx + dfyy) + f_v

            call apply_FDstencil(v_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, &
               combined_stencils, v_block%extended_block_dims, indices, v_matrix, v_rhs(global_index))

         end do
      end do

      !$omp end parallel do

   end subroutine compute_rhs

   ! Calculate velocity divergence
   subroutine calculate_velocity_divergence(u_block, v_block, u_matrix, v_matrix, FDstencil, u_x_v_y)
      type(block_type), intent(in) :: u_block, v_block
      real, dimension(:), intent(in) :: u_matrix, v_matrix
      type(FDstencil_type), target, intent(inout) :: FDstencil
      real, dimension(:), intent(inout) :: u_x_v_y

      integer :: ii, jj, global_index
      integer, dimension(u_block%ndims) :: indices, alphas
      real, dimension(:), pointer :: coefficients, dfx, dfy
      real :: u_x, v_y

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(u_block, v_block, u_matrix, v_matrix, FDstencil, u_x_v_y) &
      !$omp private(ii, jj, indices, global_index, alphas, coefficients, dfx, dfy, u_x, v_y)
      do ii = u_block%extended_block_begin_c(1)+2, u_block%extended_block_end_c(1)-1
         do jj = u_block%extended_block_begin_c(2)+2, u_block%extended_block_end_c(2)-1
            indices = [ii,jj]
            call IDX_XD(u_block%ndims, u_block%extended_block_dims, indices, global_index)

            call get_FD_coefficients_from_index(u_block%ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               u_block%extended_block_begin_c+1, u_block%extended_block_dims, indices, &
               FDstencil%scaled_stencil_coefficients, alphas, coefficients)

            dfx => coefficients(1:FDstencil%num_stencil_elements)
            dfy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call apply_FDstencil(u_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, dfx, u_block%extended_block_dims, &
               indices, u_matrix, u_x)
            call apply_FDstencil(u_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, dfy, u_block%extended_block_dims, &
               indices, v_matrix, v_y)

            u_x_v_y(global_index) = u_x + v_y

         end do
      end do

      !$omp end parallel do

   end subroutine calculate_velocity_divergence

   !> Solve the Poisson problem using the fractional step method. Currently parallel only supports Jacobi iteration.
   subroutine solve_poisson_problem(u_block, v_block, p_block, p_matrix_ptr, p_temp_matrix_ptr, u_x_v_y, &
      FDstencil, solver_params, comm, norm_array, rho, dt, converged, it)
      type(block_type), intent(inout) :: u_block, v_block, p_block
      real, dimension(:), intent(inout), pointer :: p_matrix_ptr, p_temp_matrix_ptr
      real, dimension(:), intent(in) :: u_x_v_y
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver_params
      type(comm_type), intent(in) :: comm
      real, dimension(:), intent(inout) :: norm_array
      real, intent(in) :: rho, dt
      integer, intent(out) :: converged, it

      integer :: ii, jj, p_global_index, uv_global_index
      integer, dimension(p_block%ndims) :: p_indices, uv_indices, alphas
      real, dimension(:), pointer :: coefficients, dfxx, dfyy
      real, dimension(product(FDstencil%stencil_sizes)) :: combined_stencils
      real :: local_norm, old_val, new_val, f_val

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < solver_params%max_iter)

         local_norm = 0.0

         !$omp parallel do collapse(2) default(none) &
         !$omp shared(u_block, p_block, rho, dt, u_x_v_y, p_matrix_ptr, p_temp_matrix_ptr, FDstencil) &
         !$omp private(ii, jj, p_indices, uv_indices, p_global_index, uv_global_index, old_val, alphas, coefficients, dfxx, dfyy, &
         !$omp f_val, combined_stencils, new_val) &
         !$omp reduction(+:local_norm)
         do ii = p_block%block_begin_c(1)+1, p_block%block_end_c(1)
            do jj = p_block%block_begin_c(2)+1, p_block%block_end_c(2)
               p_indices = [ii,jj]
               uv_indices = p_indices - p_block%ghost_begin
               call IDX_XD(p_block%ndims, p_block%extended_block_dims, p_indices, p_global_index)
               call IDX_XD(u_block%ndims, u_block%extended_block_dims, uv_indices, uv_global_index)

               old_val = p_matrix_ptr(p_global_index)

               call get_FD_coefficients_from_index(p_block%ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
                  p_block%extended_block_begin_c+1, p_block%extended_block_dims, p_indices, &
                  FDstencil%scaled_stencil_coefficients, alphas, coefficients)

               dfxx => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
               dfyy => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

               f_val = (rho/dt) * u_x_v_y(uv_global_index)
               combined_stencils = dfxx + dfyy

               call update_value_from_stencil(p_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, combined_stencils, &
                  p_block%extended_block_begin_c+1, p_block%extended_block_dims, p_matrix_ptr, f_val, new_val)

               p_temp_matrix_ptr(p_global_index) = new_val

               local_norm = local_norm + (abs(new_val - old_val))**2

            end do
         end do

         !$omp end parallel do

         norm_array(1) = local_norm

         call write_pressure_BC(p_block%ndims, p_block%extended_block_dims, p_block%global_begin_c, p_temp_matrix_ptr)

         call check_convergence(comm%comm, solver_params%tol, solver_params%divergence_tol, &
            1.0/product(p_block%total_grid_size), norm_array, converged)
         if(converged == -1) then
            !exit
         end if

         call swap_pointers(p_matrix_ptr, p_temp_matrix_ptr)

         call sendrecv_data_neighbors(comm%comm, p_block, p_matrix_ptr)

         it = it + 1
         converged = 0

      end do

   end subroutine solve_poisson_problem

   ! Update the velocity fields
   subroutine update_velocity_fields(u_block, v_block, p_block, FDstencil, u_star, v_star, rho, dt)
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      real, dimension(:), intent(in) :: u_star, v_star
      real, intent(in) :: rho, dt

      integer :: ii, jj, p_global_index, uv_global_index
      integer, dimension(p_block%ndims) :: p_indices, uv_indices, alphas
      real, dimension(:), pointer :: coefficients, dfx, dfy
      real, dimension(product(FDstencil%stencil_sizes)) :: combined_stencils
      real :: p_x, p_y

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(u_block, v_block, p_block, FDstencil, u_star, v_star, rho, dt) &
      !$omp private(ii, jj, p_indices, uv_indices, p_global_index, uv_global_index, alphas, coefficients, dfx, dfy, &
      !$omp combined_stencils, p_x, p_y)
      do ii = p_block%block_begin_c(1)+1, p_block%block_end_c(1)
         do jj = p_block%block_begin_c(2)+1, p_block%block_end_c(2)
            p_indices = [ii,jj]
            uv_indices = p_indices - p_block%ghost_begin
            call IDX_XD(p_block%ndims, p_block%extended_global_dims, p_indices, p_global_index)
            call IDX_XD(u_block%ndims, u_block%extended_global_dims, uv_indices, uv_global_index)

            call get_FD_coefficients_from_index(p_block%ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               p_block%extended_block_begin_c+1, p_block%extended_block_dims, p_indices, &
               FDstencil%scaled_stencil_coefficients, alphas, coefficients)

            dfx => coefficients(1:FDstencil%num_stencil_elements)
            dfy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            combined_stencils = dfx
            call apply_FDstencil(p_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, combined_stencils, &
               p_block%extended_local_size, p_indices, p_block%matrix_ptr, p_x)

            combined_stencils = dfy
            call apply_FDstencil(p_block%ndims, 1, 0, FDstencil%stencil_sizes, alphas, combined_stencils, &
               p_block%extended_local_size, p_indices, p_block%matrix_ptr, p_y)

            !write(*,*) p_indices, "P_x: ", p_x, " P_y: ", p_y

            u_block%matrix_ptr(uv_global_index) = u_star(uv_global_index) - (dt/rho) * p_x
            v_block%matrix_ptr(uv_global_index) = v_star(uv_global_index) - (dt/rho) * p_y

         end do
      end do

      !$omp end parallel do

   end subroutine update_velocity_fields

   subroutine test_pressure_bc(ndims, dims, begin_indices, p_matrix)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: dims, begin_indices
      real, dimension(:), intent(inout) :: p_matrix

      integer :: global_index
      integer, dimension(ndims) :: indices

      !$omp parallel do default(none) &
      !$omp shared(ndims, dims, begin_indices, p_matrix) &
      !$omp private(global_index, indices)
      do global_index = 1, product(dims)
         call IDX_XD_INV(ndims, dims, global_index, indices)
         p_matrix(global_index) = (indices(1)-1)*10 + indices(2)
      end do

      !$omp end parallel do

   end subroutine test_pressure_bc

end module Navier_Stokes_2D_module
