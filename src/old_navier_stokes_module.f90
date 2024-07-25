module Navier_Stokes_2D_module_old
   use constants_module, only: pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, calculate_scaled_coefficients, get_coefficients_wrapper, &
      apply_FDstencil_2D, update_value_from_stencil_2D, set_2D_matrix_coefficients
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type
   use functions_module, only: FunctionPtrType, FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers, &
      reshape_real_1D_to_2D, calculate_dt_from_CFL

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
   integer, dimension(ndims), parameter :: uv_ghost_begin = [0,0], uv_ghost_end = [0,0]
   integer, dimension(ndims), parameter :: p_ghost_begin = [0,0], p_ghost_end = [0,0]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   ! Solver parameters
   integer :: solver_type = 1, omega = 1.0
   integer, parameter :: N_iterations = 500
   integer, parameter :: N_Pressure_Poisson_iterations = 50
   real, parameter :: Pressure_Poisson_tol = 1e-6, Pressure_Poisson_div_tol = 1e-1

   ! Physical parameters
   real, parameter :: rho = 1.0, nu = 0.1, system_u_lid = 1.0, dt = 0.001
   real, dimension(ndims), parameter :: F = [0.0, 0.0]

   public :: Navier_Stokes_2D_main

contains

   !> Solve the 2D Navier-Stokes equation
   subroutine Navier_Stokes_2D_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      type(SolverParamsType):: solver_params
      type(comm_type) :: comm_params
      type(FDstencil_type) :: FDstencil_params
      type(block_type) :: u_block_params, v_block_params, p_block_params
      type(ResultType) :: result

      integer :: ii, iounit, converged, iter
      real, dimension(4) :: result_array_with_timings

      call set_SolverParamsType(Pressure_Poisson_tol, Pressure_Poisson_div_tol, N_Pressure_Poisson_iterations, 2, solver_params)

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, u_block_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, v_block_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         p_ghost_begin, p_ghost_end, stencil_begin, stencil_end, solver_type, p_block_params)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      if(solver_type == 1) then
         ! Run the solver
         do ii = 1, N_iterations
            call NS_2D_one_timestep(comm_params, u_block_params, v_block_params, p_block_params, &
               FDstencil_params, solver_params, result)
            write(*,*) "Iteration: ", ii, "/", N_iterations
            call print_resultType(result)

            if(converged == -1) then
               exit
            end if

         end do
      else
         !> Assemble the matrix A
         call assemble_matrix(p_block_params, FDstencil_params)

         ! Apply the matrix boundary conditions before decomposing the matrix
         call write_matrix_bc(p_block_params, FDstencil_params)

         !> Decompose the matrix A into LU
         call LU_decomposition(p_block_params)

         !> Write the initial condition to the system which is just zero
         !call write_inital_condition(p_block_params)

         ! do ii = 1, t_steps
         !    call NS_2D_one_timestep(comm_params, u_block_params, v_block_params, p_block_params, &
         !       FDstencil_params, solver_params, result)
         !    print*, "Time step: ", ii, " Time: ", current_time
         ! end do
      end if

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

      !call print_cart_comm_type(comm_params, iounit)

      call print_block_type(u_block_params, iounit)
      call print_block_type(v_block_params, iounit)
      call print_block_type(p_block_params, iounit)

      !call print_finite_difference_stencil(FDstencil_params, iounit)

      call close_txt_file(iounit)

      ! ! Write out system solution to a file
      if(solver_type == 1) then
         call write_block_data_to_file(u_block_params%data_layout, "output/u_solution.dat", comm_params%comm, u_block_params%matrix)
         call write_block_data_to_file(v_block_params%data_layout, "output/v_solution.dat", comm_params%comm, v_block_params%matrix)
         call write_block_data_to_file(p_block_params%data_layout, "output/p_solution.dat", comm_params%comm, p_block_params%matrix)
      else
         call write_block_data_to_file(u_block_params%data_layout, "output/u_solution.dat", comm_params%comm, u_block_params%matrix)
         call write_block_data_to_file(v_block_params%data_layout, "output/v_solution.dat", comm_params%comm, v_block_params%matrix)
         call write_block_data_to_file(p_block_params%data_layout, "output/p_solution.dat", comm_params%comm, p_block_params%matrix)
      end if

      ! Deallocate data

      call deallocate_cart_comm_type(comm_params)

      call deallocate_block_type(u_block_params)

      call deallocate_block_type(v_block_params)

      call deallocate_block_type(p_block_params)

      call deallocate_finite_difference_stencil(FDstencil_params)

   end subroutine Navier_Stokes_2D_main

   ! Solve the Navier-Stokes in 2D using the fractional step method. Ensure everything uses  [jj, ii] indexing.
   subroutine NS_2D_one_timestep(comm, u_block, v_block, p_block, FDstencil, solver, result)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver
      type(ResultType), intent(out) :: result

      ! Calculate intermediate velocity fields
      call calculate_intermediate_velocity(u_block, v_block, FDstencil)

      !> Write out the boundary conditions on the velocity fields
      call write_velocity_BC(u_block%extended_block_begin_c, u_block%extended_block_end_c, &
         u_block%f_matrix_ptr_2D, v_block%f_matrix_ptr_2D)

      ! Calculate velocity divergence
      call calculate_velocity_divergence(u_block, v_block, p_block, FDstencil)

      ! Calculate pressure correction
      call solve_poisson_problem(comm, p_block, FDstencil, solver, result)

      ! Update the velocity fields and correct the pressure
      call update_velocity_fields(u_block, v_block, p_block, FDstencil)

      ! Write out the boundary conditions on the velocity fields
      call write_velocity_BC(u_block%extended_block_begin_c, u_block%extended_block_end_c, &
         u_block%matrix_ptr_2D, v_block%matrix_ptr_2D)

   end subroutine NS_2D_one_timestep

   ! Write the lid cavity boundary conditions on the velocity fields.
   subroutine write_velocity_BC(extended_begin_c, extended_end_c, u_matrix, v_matrix)
      integer, dimension(:), intent(in) :: extended_begin_c, extended_end_c
      real, dimension(:,:), intent(inout) :: u_matrix, v_matrix

      integer :: ii, jj

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(extended_begin_c, extended_end_c, u_matrix, v_matrix) &
      !$omp private(ii, jj)
      do ii = extended_begin_c(2)+1, extended_end_c(2)
         do jj = extended_begin_c(1)+1, extended_end_c(1)
            if(jj == 1 .or. jj == extended_end_c(1) .or. ii == 1 .or. ii == extended_end_c(2)) then
               u_matrix(jj,ii) = 0.0
               v_matrix(jj,ii) = 0.0
            end if
            if(ii == 1) then
               v_matrix(jj,ii) = system_u_lid
            end if
         end do
      end do

      !$omp end parallel do

   end subroutine write_velocity_BC

   ! Write the lid cavity boundary conditions on the pressure field.
   subroutine write_pressure_BC(extended_begin_c, extended_end_c, p_matrix)
      integer, dimension(:), intent(in) :: extended_begin_c, extended_end_c
      real, dimension(:,:), intent(inout) :: p_matrix

      integer :: ii, jj

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(extended_begin_c, extended_end_c, p_matrix) &
      !$omp private(ii, jj)
      do ii = extended_begin_c(2)+1, extended_end_c(2)
         do jj = extended_begin_c(1)+1, extended_end_c(1)

            ! Left wall: u_x = 0 (Neumann condition)
            if(jj == 1) then
               p_matrix(jj,ii) = p_matrix(jj+1,ii)
            end if

            ! Right wall: u_x = 0 (Neumann condition)
            if(jj == extended_end_c(1)) then
               p_matrix(jj,ii) = p_matrix(jj-1,ii)
            end if

            ! Top wall: p = 0 (Dirichlet condition)
            if (ii == 1) then
               p_matrix(jj,ii) = 0.0
            end if

            ! Bottom wall: u_y = 0 (Neumann condition)
            if(ii == extended_end_c(2)) then
               p_matrix(jj,ii) = p_matrix(jj,ii-1)
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_pressure_BC

   ! Calculate rhs for the 2D Navier-Stokes equation
   subroutine calculate_intermediate_velocity(u_block, v_block, FDstencil)
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
      !$omp shared(u_block, v_block, FDstencil) &
      !$omp private(ii, jj, uv_local_indices, u, v, u_rhs, v_rhs, coefficients, alpha, beta, dx, dy, dxx, dyy, combined_stencil, &
      !$omp dx_2D, dy_2D, dxx_2D, dyy_2D, combined_stencil_2D)
      do ii = u_block%block_begin_c(2)+2, u_block%block_end_c(2)-1
         do jj = u_block%block_begin_c(1)+2, u_block%block_end_c(1)-1
            uv_local_indices = [jj,ii]

            u = u_block%matrix_ptr_2D(jj,ii)
            v = v_block%matrix_ptr_2D(jj,ii)

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
            u_block%f_matrix_ptr_2D(jj,ii) = u_block%matrix_ptr_2D(jj,ii) + dt * u_rhs
            v_block%f_matrix_ptr_2D(jj,ii) = v_block%matrix_ptr_2D(jj,ii) + dt * v_rhs

         end do
      end do

      !$omp end parallel do

   end subroutine calculate_intermediate_velocity

   ! Calculate velocity divergence
   subroutine calculate_velocity_divergence(u_block, v_block, p_block, FDstencil)
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(ndims) :: uv_local_indices, p_local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dx, dy
      real, contiguous, dimension(:,:), pointer :: dx_2D, dy_2D
      real :: u_x, v_y

      call calculate_scaled_coefficients(u_block%ndims, u_block%extended_grid_dx, FDstencil)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(u_block, v_block, p_block, FDstencil) &
      !$omp private(ii, jj, uv_local_indices, p_local_indices, alpha, beta, coefficients, dx, dy, dx_2D, dy_2D, u_x, v_y)
      do ii = u_block%block_begin_c(2)+2, u_block%block_end_c(2)-1
         do jj = u_block%block_begin_c(1)+2, u_block%block_end_c(1)-1
            uv_local_indices = [jj,ii]
            p_local_indices = p_block%block_begin_c + uv_local_indices - u_block%block_begin_c

            call get_coefficients_wrapper(FDstencil, [1,1], u_block%extended_block_dims, &
               uv_local_indices, alpha, beta, coefficients)

            dy => coefficients(1:FDstencil%num_stencil_elements)
            dx => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dy, dy_2D)

            call apply_FDstencil_2D(dx_2D, u_block%f_matrix_ptr_2D, uv_local_indices, alpha, beta, u_x)
            call apply_FDstencil_2D(dy_2D, v_block%f_matrix_ptr_2D, uv_local_indices, alpha, beta, v_y)

            p_block%f_matrix_ptr_2D(p_local_indices(1), p_local_indices(2)) = (rho/dt) * (u_x + v_y)

         end do
      end do

      !$omp end parallel do

   end subroutine calculate_velocity_divergence

   !> Solve the Poisson problem for the pressure field
   subroutine solve_poisson_problem(comm, p_block, FDstencil, solver, result)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver
      type(ResultType), intent(out) :: result

      integer :: ii, jj, it, converged
      integer, dimension(ndims) :: p_index, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencil
      real, contiguous, dimension(:,:), pointer :: combined_stencil_2D
      real :: u0_val, u1_val, f_val, r1_val, local_norm

      call calculate_scaled_coefficients(p_block%ndims, p_block%extended_grid_dx, FDstencil)
      call write_pressure_BC(p_block%extended_block_begin_c, p_block%extended_block_end_c, p_block%matrix_ptr_2D)
      ! SendRECV if needed

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < solver%max_iter)

         local_norm = 0.0
         do ii = p_block%block_begin_c(2)+2, p_block%block_end_c(2)-1
            do jj = p_block%block_begin_c(1)+2, p_block%block_end_c(1)-1
               p_index = [ii,jj]

               f_val = p_block%f_matrix_ptr_2D(jj,ii)
               u0_val = p_block%matrix_ptr_2D(jj,ii)

               call get_coefficients_wrapper(FDstencil, [1,1], p_block%extended_block_end_c, &
                  p_index, alpha, beta, coefficients)

               dyy => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
               dxx => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

               combined_stencil = dxx + dyy

               call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencil, combined_stencil_2D)

               call update_value_from_stencil_2D(combined_stencil_2D, p_block%matrix_ptr_2D, p_index, alpha, beta, &
                  f_val, u1_val, r1_val)

               u1_val = (1.0 - omega) * u0_val + omega * u1_val

               p_block%matrix_ptr_2D(jj,ii) = u1_val

               local_norm = local_norm + (abs(u1_val - u0_val))**2

            end do
         end do

         !norm_array(1) = local_norm

         call write_pressure_BC(p_block%extended_block_begin_c, p_block%extended_block_end_c, p_block%matrix_ptr_2D)
         ! SendRECV if needed

         !call check_convergence(comm, tol, 1e6, 1.0/product(p_dims), norm_array, converged)
         if(converged == -1) then
            !exit
         end if

         !call swap_pointers(p_ptr, p_temp_ptr)

         it = it + 1

      end do

      call set_ResultType(converged, it, local_norm, local_norm, local_norm, local_norm, result)

   end subroutine solve_poisson_problem

   ! Update the velocity fields using the pressure correction
   subroutine update_velocity_fields(u_block, v_block, p_block, FDstencil)
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(ndims) :: uv_local_indices, p_local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dx, dy
      real, contiguous, dimension(:,:), pointer :: dx_2D, dy_2D
      real :: p_x, p_y

      call calculate_scaled_coefficients(p_block%ndims, p_block%extended_grid_dx, FDstencil)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(u_block, v_block, p_block, FDstencil) &
      !$omp private(ii, jj, uv_local_indices, p_local_indices, alpha, beta, coefficients, dx, dy, dx_2D, dy_2D, p_x, p_y)
      do ii = u_block%block_begin_c(2)+2, u_block%block_end_c(2)-1
         do jj = u_block%block_begin_c(1)+2, u_block%block_end_c(1)-1
            uv_local_indices = [jj,ii]
            p_local_indices = p_block%block_begin_c + uv_local_indices - u_block%block_begin_c

            call get_coefficients_wrapper(FDstencil, [1,1], p_block%extended_block_dims, &
               p_local_indices, alpha, beta, coefficients)

            dy => coefficients(1:FDstencil%num_stencil_elements)
            dx => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dy, dy_2D)

            call apply_FDstencil_2D(dx_2D, p_block%matrix_ptr_2D, p_local_indices, alpha, beta, p_x)
            call apply_FDstencil_2D(dy_2D, p_block%matrix_ptr_2D, p_local_indices, alpha, beta, p_y)

            u_block%matrix_ptr_2D(jj,ii) = u_block%f_matrix_ptr_2D(jj,ii) - (dt/rho) * p_x
            v_block%matrix_ptr_2D(jj,ii) = v_block%f_matrix_ptr_2D(jj,ii) - (dt/rho) * p_y

         end do
      end do

      !$omp end parallel do

   end subroutine update_velocity_fields

   !> Assembly of matrix that represents the 2D Navier-Stokes equation
   subroutine assemble_matrix(block_params, FDstencil_params)
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), intent(inout) :: FDstencil_params

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil_params%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D

      integer :: num_equations, num_rhs, lda, ldb
      integer, dimension(:), allocatable :: ipiv
      integer :: info

      call calculate_scaled_coefficients(block_params%ndims, block_params%extended_grid_dx, FDstencil_params)

      do ii = block_params%extended_block_begin_c(2)+1, block_params%extended_block_end_c(2)
         do jj = block_params%extended_block_begin_c(1)+1, block_params%extended_block_end_c(1)
            local_indices = [jj,ii]

            call get_coefficients_wrapper(FDstencil_params, [1,1], block_params%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(1:FDstencil_params%num_stencil_elements)
            dyy => coefficients(FDstencil_params%num_stencil_elements + 1:2 * FDstencil_params%num_stencil_elements)

            combined_stencils = dxx + dyy

            call reshape_real_1D_to_2D(FDstencil_params%stencil_sizes, combined_stencils, combined_stencils_2D)

            call set_2D_matrix_coefficients(block_params%extended_block_dims, combined_stencils_2D, &
               block_params%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)

         end do
      end do

   end subroutine assemble_matrix

   !> Write the Dirichlet and Neumann boundary conditions to the matrix
   subroutine write_matrix_bc(p_block, FDstencil)
      type(block_type), intent(inout) :: p_block
      type(FDstencil_type), intent(inout) :: FDstencil

      integer :: ii, jj, diag
      integer, dimension(ndims) :: local_indices, global_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients
      real, contiguous, dimension(:), pointer :: dx, dy
      real, contiguous, dimension(:,:), pointer :: dx_2D, dy_2D

      do ii = p_block%extended_block_begin_c(2)+1, p_block%extended_block_end_c(2)
         do jj = p_block%extended_block_begin_c(1)+1, p_block%extended_block_end_c(1)
            local_indices = [jj,ii]
            global_indices = p_block%extended_global_begin_c + local_indices

            diag = local_indices(1) + (local_indices(2) - 1) * p_block%extended_block_dims(1)

            call get_coefficients_wrapper(FDstencil, [1,1], p_block%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dy => coefficients(1:FDstencil%num_stencil_elements)
            dx => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dy, dy_2D)

            ! At the top wall, p = 0 (Dirichlet condition)
            if(ii == 1) then
               p_block%direct_solver_matrix_ptr_2D(diag,:) = 0.0
               p_block%direct_solver_matrix_ptr_2D(diag,diag) = 1.0
            end if

            ! At the bottom wall, u_y = 0 (Neumann condition)
            if(ii == p_block%extended_block_end_c(2)) then
               call set_2D_matrix_coefficients(p_block%extended_block_dims, dy_2D, &
                  p_block%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)
            end if

            ! At the left and right walls, u_x = 0 (Neumann condition)
            if(jj == 1 .or. jj == p_block%extended_block_end_c(1)) then
               call set_2D_matrix_coefficients(p_block%extended_block_dims, dx_2D, &
                  p_block%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)
            end if

         end do
      end do

   end subroutine write_matrix_bc

end module Navier_Stokes_2D_module_old
