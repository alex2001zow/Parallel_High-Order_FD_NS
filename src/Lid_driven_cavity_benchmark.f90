module Lid_driven_cavity_benchmark_module
   use constants_module, only: pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, calculate_scaled_coefficients, get_coefficients_wrapper, &
      apply_FDstencil_2D, update_value_from_stencil_2D, set_matrix_coefficients
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type
   use functions_module, only: FunctionPtrType, FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers, &
      reshape_real_1D_to_2D, calculate_dt_from_CFL, swap_pointers_2D, copy_1D_array

   use solver_module, only: SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, check_convergence, set_ResultType, LU_decomposition, solve_LU_system
   implicit none

   private

   !> Number of dimensions and number of derivatives
   integer, parameter :: ndims = 2
   integer, parameter :: num_derivatives = 4
   integer, dimension(ndims*num_derivatives), parameter :: derivatives = [1,0,0,1,2,0,0,2] ! dy, dx, dyy, dxx

   ! Grid parameters
   integer, dimension(ndims) :: grid_size = [42,42], processor_dims = [1,1]
   logical, dimension(ndims) :: periods = [.false.,.false.]
   logical, parameter :: reorder = .true.
   real, dimension(ndims) :: domain_begin = [0,0], domain_end = [1,1]
   integer, dimension(ndims), parameter :: stencil_sizes = 3
   integer, dimension(ndims), parameter :: uv_ghost_begin = [1,1], uv_ghost_end = [1,1]
   integer, dimension(ndims), parameter :: p_ghost_begin = [1,1], p_ghost_end = [1,1]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   ! Solver parameters
   integer, parameter :: direct_or_iterative = 1, Jacobi_or_GS = 1, omega = 1.0
   integer, parameter :: N_iterations = 500
   integer, parameter :: N_Pressure_Poisson_iterations = 50
   real, parameter :: Pressure_Poisson_tol = 1e-6, Pressure_Poisson_div_tol = 1e-1

   ! Physical parameters
   real, parameter :: rho = 1.0, nu = 0.1, system_u_lid = 1.0
   real, dimension(ndims), parameter :: F = [0.0, 0.0]

   public :: Lid_driven_cavity_benchmark_2D

contains

   !> Solve the 2D Navier-Stokes equation
   subroutine Lid_driven_cavity_benchmark_2D(rank, world_size)
      integer, intent(in) :: rank, world_size

      type(SolverParamsType):: solver
      type(comm_type) :: comm
      type(FDstencil_type) :: FDstencil
      type(block_type) :: u1_block, v1_block, p1_block, &
         u2_block, v2_block, p2_block, &
         u3_block, v3_block, p3_block, &
         u4_block, v4_block, p4_block
      type(ResultType) :: result

      integer :: ii, iounit
      real, dimension(4) :: result_array_with_timings
      real :: dt

      call set_SolverParamsType(Pressure_Poisson_tol, Pressure_Poisson_div_tol, N_Pressure_Poisson_iterations, Jacobi_or_GS, &
         solver)

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, u1_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, v1_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         p_ghost_begin, p_ghost_end, stencil_begin, stencil_end, direct_or_iterative, p1_block)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, u2_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, v2_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         p_ghost_begin, p_ghost_end, stencil_begin, stencil_end, direct_or_iterative, p2_block)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, u3_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, v3_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         p_ghost_begin, p_ghost_end, stencil_begin, stencil_end, direct_or_iterative, p3_block)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, u4_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, v4_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         p_ghost_begin, p_ghost_end, stencil_begin, stencil_end, direct_or_iterative, p4_block)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil)

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      dt = 0.001

      if(direct_or_iterative == 1) then
         ! Run the solver
         do ii = 1, N_iterations
            call RK4(dt, comm, u1_block, v1_block, p1_block, &
               u2_block, v2_block, p2_block, &
               u3_block, v3_block, p3_block, &
               u4_block, v4_block, p4_block, &
               FDstencil, solver, result)
            if(rank == MASTER_RANK) then
               print*, "Iteration: ", ii, "/", N_iterations
               call print_resultType(result)
            end if
            if(result%converged == -1) then
               exit
            end if

         end do
      else
         !> Assemble the matrix A
         call assemble_matrix_2D(p1_block, FDstencil)

         ! Apply the matrix boundary conditions before decomposing the matrix
         call write_matrix_bc_2D(p1_block, FDstencil)

         !> Decompose the matrix A into LU
         call LU_decomposition(p1_block)

         !> Write the initial condition to the system which is just zero

         do ii = 1, N_iterations
            call one_timestep_2D(dt, 1.0, comm, u1_block, v1_block, p1_block, FDstencil, solver, result)
            if(rank == MASTER_RANK) then
               print*, "Iteration: ", ii, "/", N_iterations
               call print_resultType(result)
            end if
         end do
      end if

      result_array_with_timings(2) = MPI_WTIME()

      result_array_with_timings(3) = result_array_with_timings(2) - result_array_with_timings(1)

      call all_reduce_mpi_wrapper(result_array_with_timings(3), result_array_with_timings(4), 1, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm%comm)

      ! Write out the result and timings from the master rank
      if(rank == MASTER_RANK) then
         write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(4)/world_size, " seconds"
      end if

      ! Write out our setup to a file. Just for debugging purposes
      !call open_txt_file("output/output_from_", rank, iounit)

      !call print_cart_comm_type(comm_params, iounit)

      !call print_block_type(u_block_params, iounit)
      !call print_block_type(v_block_params, iounit)
      !call print_block_type(p_block_params, iounit)

      !call print_finite_difference_stencil(FDstencil_params, iounit)

      !call close_txt_file(iounit)

      ! ! Write out system solution to a file
      if(direct_or_iterative == 1) then
         call write_block_data_to_file(u1_block%data_layout, "output/u_solution.dat", comm%comm, u1_block%matrix)
         call write_block_data_to_file(v1_block%data_layout, "output/v_solution.dat", comm%comm, v1_block%matrix)
         call write_block_data_to_file(p1_block%data_layout, "output/p_solution.dat", comm%comm, p1_block%matrix)
      else
         call write_block_data_to_file(u1_block%data_layout, "output/u_solution.dat", comm%comm, u1_block%matrix)
         call write_block_data_to_file(v1_block%data_layout, "output/v_solution.dat", comm%comm, v1_block%matrix)
         call write_block_data_to_file(p1_block%data_layout, "output/p_solution.dat", comm%comm, p1_block%matrix)
      end if

      ! Deallocate data

      call deallocate_cart_comm_type(comm)

      call deallocate_block_type(u1_block)
      call deallocate_block_type(v1_block)
      call deallocate_block_type(p1_block)

      call deallocate_block_type(u2_block)
      call deallocate_block_type(v2_block)
      call deallocate_block_type(p2_block)

      call deallocate_block_type(u3_block)
      call deallocate_block_type(v3_block)
      call deallocate_block_type(p3_block)

      call deallocate_block_type(u4_block)
      call deallocate_block_type(v4_block)
      call deallocate_block_type(p4_block)

      call deallocate_finite_difference_stencil(FDstencil)

   end subroutine Lid_driven_cavity_benchmark_2D

   !> Subroutine to timestep using the RK4-method
   subroutine RK4(dt, comm, u1_block, v1_block, p1_block, &
      u2_block, v2_block, p2_block, &
      u3_block, v3_block, p3_block, &
      u4_block, v4_block, p4_block, &
      FDstencil, solver, result)
      real, intent(in) :: dt
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: u1_block, v1_block, p1_block, &
         u2_block, v2_block, p2_block, &
         u3_block, v3_block, p3_block, &
         u4_block, v4_block, p4_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver
      type(ResultType), intent(out) :: result

      !> RK4 time coefficients
      real, parameter :: RK1_time = 0, RK2_time = 1.0/2.0, RK3_time = 1.0/2.0, RK4_time = 1.0

      !> RK4 combination coefficients
      real, parameter :: RK1_coef = 1.0/6.0, RK2_coef = 1.0/3.0, RK3_coef = 1.0/3.0, RK4_coef = 1.0/6.0

      !> RK4 step 1 is just the current state
      ! Do nothing

      !> RK4 step 2
      call one_timestep_2D(dt, RK2_time, comm, u2_block, v2_block, p2_block, FDstencil, solver, result)

      !> RK4 step 3
      call one_timestep_2D(dt, RK3_time, comm, u3_block, v3_block, p3_block, FDstencil, solver, result)

      !> RK4 step 4
      call one_timestep_2D(dt, RK4_time, comm, u4_block, v4_block, p4_block, FDstencil, solver, result)

      !> Combine the results
      u1_block%matrix_ptr = u1_block%matrix_ptr*RK1_coef + u2_block%matrix_ptr*RK2_coef + &
         u3_block%matrix_ptr*RK3_coef + u4_block%matrix_ptr*RK4_coef

      v1_block%matrix_ptr = v1_block%matrix_ptr*RK1_coef + v2_block%matrix_ptr*RK2_coef + &
         v3_block%matrix_ptr*RK3_coef + v4_block%matrix_ptr*RK4_coef

      p1_block%matrix_ptr = p1_block%matrix_ptr*RK1_coef + p2_block%matrix_ptr*RK2_coef + &
         p3_block%matrix_ptr*RK3_coef + p4_block%matrix_ptr*RK4_coef

   end subroutine RK4

   ! Solve the Navier-Stokes in 2D using the fractional step method.
   subroutine one_timestep_2D(dt, RK_coef, comm, u_block, v_block, p_block, FDstencil, solver, result)
      real, intent(in) :: dt, RK_coef
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver
      type(ResultType), intent(out) :: result

      ! Calculate intermediate velocity fields
      call calculate_intermediate_velocity_2D(dt, RK_coef, u_block, v_block, FDstencil)

      !> Write out the boundary conditions on the velocity fields
      call write_velocity_BC_2D(u_block%extended_block_begin_c, u_block%extended_block_end_c, &
         u_block%f_matrix_ptr_2D, v_block%f_matrix_ptr_2D)

      ! Calculate velocity divergence
      call calculate_velocity_divergence_2D(dt, u_block, v_block, p_block, FDstencil)

      if(direct_or_iterative == 1) then
         ! Calculate pressure correction
         call solve_poisson_problem_2D(comm, p_block, FDstencil, solver, result)
      else
         p_block%f_matrix_ptr_2D(:,1) = 0.0
         p_block%f_matrix_ptr_2D(:,p_block%extended_block_end_c(2)) = 0.0
         p_block%f_matrix_ptr_2D(1,:) = 0.0
         p_block%f_matrix_ptr_2D(p_block%extended_block_end_c(1),:) = 0.0
         call solve_LU_system(p_block)
         call copy_1D_array(p_block%f_matrix_ptr, p_block%matrix_ptr)
      end if

      ! Update the velocity fields and correct the pressure
      call update_velocity_fields_2D(dt, u_block, v_block, p_block, FDstencil)

      ! Write out the boundary conditions on the velocity fields
      call write_velocity_BC_2D(u_block%extended_block_begin_c, u_block%extended_block_end_c, &
         u_block%matrix_ptr_2D, v_block%matrix_ptr_2D)

   end subroutine one_timestep_2D

   ! Write the lid cavity boundary conditions on the velocity fields.
   subroutine write_velocity_BC_2D(extended_begin_c, extended_end_c, u_matrix, v_matrix)
      integer, dimension(:), intent(in) :: extended_begin_c, extended_end_c
      real, dimension(:,:), intent(inout) :: u_matrix, v_matrix

      integer :: ii, jj
      integer, dimension(2) :: local_indices

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(extended_begin_c, extended_end_c, u_matrix, v_matrix) &
      !$omp private(ii, jj, local_indices)
      do jj = extended_begin_c(2)+1, extended_end_c(2)
         do ii = extended_begin_c(1)+1, extended_end_c(1)
            local_indices = [ii,jj]

            ! All boundaries are no-slip
            if(local_indices(1) == 1 .or. local_indices(1) == extended_end_c(1) .or. &
               local_indices(2) == 1 .or. local_indices(2) == extended_end_c(2)) then
               u_matrix(local_indices(1),local_indices(2)) = 0.0
               v_matrix(local_indices(1),local_indices(2)) = 0.0
            end if

            ! Except the top wall which is moving
            if(local_indices(2) == 1) then
               u_matrix(local_indices(1),local_indices(2)) = system_u_lid
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_velocity_BC_2D

   ! Write the lid cavity boundary conditions on the pressure field.
   subroutine write_pressure_BC_2D(extended_begin_c, extended_end_c, p_matrix)
      integer, dimension(:), intent(in) :: extended_begin_c, extended_end_c
      real, dimension(:,:), intent(inout) :: p_matrix

      integer :: ii, jj
      integer, dimension(2) :: local_indices

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(extended_begin_c, extended_end_c, p_matrix) &
      !$omp private(ii, jj, local_indices)
      do jj = extended_begin_c(2)+1, extended_end_c(2)
         do ii = extended_begin_c(1)+1, extended_end_c(1)
            local_indices = [ii,jj]

            ! Left wall: u_x = 0 (Neumann condition)
            if(local_indices(1) == 1) then
               p_matrix(local_indices(1),local_indices(2)) = p_matrix(local_indices(1)+2,local_indices(2))
            end if

            ! Right wall: u_x = 0 (Neumann condition)
            if(local_indices(1) == extended_end_c(1)) then
               p_matrix(local_indices(1),local_indices(2)) = p_matrix(local_indices(1)-2,local_indices(2))
            end if

            ! Top wall: p = 0 (Dirichlet condition)
            if (local_indices(2) == 1) then
               p_matrix(local_indices(1),local_indices(2)) = 0.0
            end if

            ! Bottom wall: u_y = 0 (Neumann condition)
            if(local_indices(2) == extended_end_c(2)) then
               p_matrix(local_indices(1),local_indices(2)) = p_matrix(local_indices(1),local_indices(2)-2)
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_pressure_BC_2D

   ! Calculate rhs for the 2D Navier-Stokes equation
   subroutine calculate_intermediate_velocity_2D(dt, RK_coef, u_block, v_block, FDstencil)
      real, intent(in) :: dt, RK_coef
      type(block_type), intent(inout) :: u_block, v_block
      type(FDstencil_type), target, intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(ndims) :: uv_local_indices, alpha, beta
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencil
      real, contiguous, dimension(:), pointer :: coefficients, dx, dy, dxx, dyy
      real, contiguous, dimension(:,:), pointer :: dx_2D, dy_2D, dxx_2D, dyy_2D, combined_stencil_2D
      real :: u, v, u_rhs, v_rhs

      call calculate_scaled_coefficients(u_block%ndims, u_block%extended_grid_dx, FDstencil)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(dt, RK_coef, u_block, v_block, FDstencil) &
      !$omp private(ii, jj, uv_local_indices, u, v, u_rhs, v_rhs, coefficients, alpha, beta, dx, dy, dxx, dyy, combined_stencil, &
      !$omp dx_2D, dy_2D, dxx_2D, dyy_2D, combined_stencil_2D)
      do jj = u_block%block_begin_c(2)+1, u_block%block_end_c(2)
         do ii = u_block%block_begin_c(1)+1, u_block%block_end_c(1)
            uv_local_indices = [ii,jj]

            u = u_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2))
            v = v_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2))

            call get_coefficients_wrapper(FDstencil, [1,1], u_block%extended_block_dims, &
               uv_local_indices, alpha, beta, coefficients)

            dy => coefficients(1:FDstencil%num_stencil_elements)
            dx => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)
            dyy => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
            dxx => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

            combined_stencil = -(u*dx + v*dy) + nu * (dxx + dyy) + F(1)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencil, combined_stencil_2D)
            call apply_FDstencil_2D(combined_stencil_2D, u_block%matrix_ptr_2D, uv_local_indices, alpha, beta, u_rhs)

            combined_stencil = -(u*dx + v*dy) + nu * (dxx + dyy) + F(2)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencil, combined_stencil_2D)
            call apply_FDstencil_2D(combined_stencil_2D, v_block%matrix_ptr_2D, uv_local_indices, alpha, beta, v_rhs)

            ! Euler time step
            u_block%f_matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2)) = &
               u_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2)) + RK_coef * dt * u_rhs

            v_block%f_matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2)) = &
               v_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2)) + RK_coef * dt * v_rhs

         end do
      end do

      !$omp end parallel do

   end subroutine calculate_intermediate_velocity_2D

   ! Calculate velocity divergence
   subroutine calculate_velocity_divergence_2D(dt, u_block, v_block, p_block, FDstencil)
      real, intent(in) :: dt
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(ndims) :: uv_local_indices, p_local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dx, dy
      real, contiguous, dimension(:,:), pointer :: dx_2D, dy_2D
      real :: u_x, v_y

      call calculate_scaled_coefficients(u_block%ndims, u_block%extended_grid_dx, FDstencil)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(dt, u_block, v_block, p_block, FDstencil) &
      !$omp private(ii, jj, uv_local_indices, p_local_indices, alpha, beta, coefficients, dx, dy, dx_2D, dy_2D, u_x, v_y)
      do jj = u_block%block_begin_c(2)+1, u_block%block_end_c(2)
         do ii = u_block%block_begin_c(1)+1, u_block%block_end_c(1)
            uv_local_indices = [ii,jj]
            p_local_indices = p_block%block_begin_c + uv_local_indices - u_block%block_begin_c

            call get_coefficients_wrapper(FDstencil, [1,1], u_block%extended_block_dims, &
               uv_local_indices, alpha, beta, coefficients)

            dy => coefficients(1:FDstencil%num_stencil_elements)
            dx => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dy, dy_2D)

            call apply_FDstencil_2D(dx_2D, v_block%f_matrix_ptr_2D, uv_local_indices, alpha, beta, u_x)
            call apply_FDstencil_2D(dy_2D, u_block%f_matrix_ptr_2D, uv_local_indices, alpha, beta, v_y)

            p_block%f_matrix_ptr_2D(p_local_indices(1), p_local_indices(2)) = (rho/dt) * (u_x + v_y)

         end do
      end do

      !$omp end parallel do

   end subroutine calculate_velocity_divergence_2D

   !> Solve the Poisson problem for the pressure field
   subroutine solve_poisson_problem_2D(comm, p_block, FDstencil, solver, result)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver
      type(ResultType), intent(out) :: result

      integer :: ii, jj, it, converged
      integer, dimension(ndims) :: p_local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencil
      real, contiguous, dimension(:,:), pointer :: combined_stencil_2D
      real :: u0_val, u1_val, f_val, r1_val, local_norm

      call calculate_scaled_coefficients(p_block%ndims, p_block%extended_grid_dx, FDstencil)
      call write_pressure_BC_2D(p_block%extended_block_begin_c, p_block%extended_block_end_c, p_block%matrix_ptr_2D)
      ! SendRECV if needed

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < solver%max_iter)

         local_norm = 0.0
         !$omp parallel do collapse(2) default(none) reduction(+:local_norm) &
         !$omp shared(p_block, FDstencil) &
         !$omp private(ii, jj, p_local_indices, alpha, beta, coefficients, dxx, dyy, combined_stencil, combined_stencil_2D, &
         !$omp u0_val, u1_val, f_val, r1_val)
         do jj = p_block%block_begin_c(2)+1, p_block%block_end_c(2)
            do ii = p_block%block_begin_c(1)+1, p_block%block_end_c(1)
               p_local_indices = [ii,jj]

               f_val = p_block%f_matrix_ptr_2D(p_local_indices(1),p_local_indices(2))
               u0_val = p_block%matrix_ptr_2D(p_local_indices(1),p_local_indices(2))

               call get_coefficients_wrapper(FDstencil, [1,1], p_block%extended_block_end_c, &
                  p_local_indices, alpha, beta, coefficients)

               dyy => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
               dxx => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

               combined_stencil = dxx + dyy
               call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencil, combined_stencil_2D)
               call update_value_from_stencil_2D(combined_stencil_2D, p_block%matrix_ptr_2D, p_local_indices, alpha, beta, &
                  f_val, u1_val, r1_val)

               u1_val = (1.0 - omega) * u0_val + omega * u1_val

               p_block%temp_matrix_ptr_2D(p_local_indices(1),p_local_indices(2)) = u1_val

               local_norm = local_norm + (abs(u1_val - u0_val))**2

            end do
         end do

         !$omp end parallel do

         call swap_pointers_2D(p_block%matrix_ptr_2D, p_block%temp_matrix_ptr_2D)

         call write_pressure_BC_2D(p_block%extended_block_begin_c, p_block%extended_block_end_c, p_block%matrix_ptr_2D)
         ! SendRECV if needed

         !norm_array(1) = local_norm

         !call check_convergence(comm, tol, 1e6, 1.0/product(p_dims), norm_array, converged)
         if(converged == -1) then
            !exit
         end if

         it = it + 1

      end do

      call set_ResultType(converged, it, local_norm, local_norm, local_norm, local_norm, result)

   end subroutine solve_poisson_problem_2D

   ! Update the velocity fields using the pressure correction
   subroutine update_velocity_fields_2D(dt, u_block, v_block, p_block, FDstencil)
      real, intent(in) :: dt
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(ndims) :: uv_local_indices, p_local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dx, dy
      real, contiguous, dimension(:,:), pointer :: dx_2D, dy_2D
      real :: p_x, p_y

      call calculate_scaled_coefficients(p_block%ndims, p_block%extended_grid_dx, FDstencil)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(dt, u_block, v_block, p_block, FDstencil) &
      !$omp private(ii, jj, uv_local_indices, p_local_indices, alpha, beta, coefficients, dx, dy, dx_2D, dy_2D, p_x, p_y)
      do jj = u_block%block_begin_c(2)+1, u_block%block_end_c(2)
         do ii = u_block%block_begin_c(1)+1, u_block%block_end_c(1)
            uv_local_indices = [ii,jj]
            p_local_indices = p_block%block_begin_c + uv_local_indices - u_block%block_begin_c

            call get_coefficients_wrapper(FDstencil, [1,1], p_block%extended_block_dims, &
               p_local_indices, alpha, beta, coefficients)

            dy => coefficients(1:FDstencil%num_stencil_elements)
            dx => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dy, dy_2D)

            call apply_FDstencil_2D(dx_2D, p_block%matrix_ptr_2D, p_local_indices, alpha, beta, p_x)
            call apply_FDstencil_2D(dy_2D, p_block%matrix_ptr_2D, p_local_indices, alpha, beta, p_y)

            v_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2)) = &
               v_block%f_matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2))  - (dt/rho) * p_x
            u_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2))  = &
               u_block%f_matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2))  - (dt/rho) * p_y

         end do
      end do

      !$omp end parallel do

   end subroutine update_velocity_fields_2D

   !> Assembly of matrix that represents the 2D Navier-Stokes equation
   subroutine assemble_matrix_2D(p_block, FDstencil)
      type(block_type), intent(inout) :: p_block
      type(FDstencil_type), intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils

      call calculate_scaled_coefficients(p_block%ndims, p_block%extended_grid_dx, FDstencil)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(p_block, FDstencil) &
      !$omp private(ii, jj, local_indices, alpha, beta, coefficients, dxx, dyy, combined_stencils)
      do jj = p_block%extended_block_begin_c(2)+1, p_block%extended_block_end_c(2)
         do ii = p_block%extended_block_begin_c(1)+1, p_block%extended_block_end_c(1)
            local_indices = [ii,jj]

            call get_coefficients_wrapper(FDstencil, [1,1], p_block%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dyy => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
            dxx => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dyy

            call set_matrix_coefficients(p_block%ndims, FDstencil%stencil_sizes, p_block%extended_block_dims, &
               combined_stencils, p_block%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)

         end do
      end do

      !$omp end parallel do

   end subroutine assemble_matrix_2D

   !> Write the Dirichlet and Neumann boundary conditions to the matrix
   subroutine write_matrix_bc_2D(p_block, FDstencil)
      type(block_type), intent(inout) :: p_block
      type(FDstencil_type), intent(inout) :: FDstencil

      integer :: ii, jj, diag
      integer, dimension(ndims) :: local_indices, global_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients
      real, contiguous, dimension(:), pointer :: dx, dy

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(p_block, FDstencil) &
      !$omp private(ii, jj, local_indices, global_indices, diag, alpha, beta, coefficients, dx, dy)
      do jj = p_block%extended_block_begin_c(2)+1, p_block%extended_block_end_c(2)
         do ii = p_block%extended_block_begin_c(1)+1, p_block%extended_block_end_c(1)
            local_indices = [ii,jj]
            global_indices = p_block%extended_global_begin_c + local_indices

            call IDX_XD(p_block%ndims, p_block%extended_block_dims, local_indices, diag)

            call get_coefficients_wrapper(FDstencil, [1,1], p_block%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dy => coefficients(1:FDstencil%num_stencil_elements)
            dx => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            ! At the left and right walls, u_x = 0 (Neumann condition)
            if(local_indices(1) == 1 .or. local_indices(1) == p_block%extended_block_end_c(1)) then
               call set_matrix_coefficients(p_block%ndims, FDstencil%stencil_sizes, p_block%extended_block_dims, &
                  dy, p_block%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)
            end if

            ! At the top wall, p = 0 (Dirichlet condition)
            if(local_indices(2) == 1) then
               p_block%direct_solver_matrix_ptr_2D(diag,:) = 0.0
               p_block%direct_solver_matrix_ptr_2D(diag,diag) = 1.0
            end if

            ! At the bottom wall, u_y = 0 (Neumann condition)
            if(local_indices(2) == p_block%extended_block_end_c(2)) then
               call set_matrix_coefficients(p_block%ndims, FDstencil%stencil_sizes, p_block%extended_block_dims, &
                  dx, p_block%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_matrix_bc_2D

end module Lid_driven_cavity_benchmark_module
