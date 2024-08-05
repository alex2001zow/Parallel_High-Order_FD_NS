module flow_around_cyliner_module
   use constants_module, only: pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, calculate_scaled_coefficients, get_coefficients_wrapper, &
      apply_FDstencil_2D, update_value_from_stencil_2D, set_matrix_coefficients
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type,sendrecv_data_neighbors
   use functions_module, only: calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers, &
      reshape_real_1D_to_2D, calculate_dt_from_CFL, swap_pointers_2D, LU_decomposition, solve_LU_system, &
      copy_vector, scale_vector, daxpy_to_vector, set_bc_zero_2D

   use solver_module, only: SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, check_convergence, set_ResultType

   use multigrid_module, only: full_weighing_restriction_2D, bilinear_prolongation_2D, apply_correction
   use Poisson_functions_module, only: Poisson_Gauss_Seidel_2D, Poisson_Gauss_Seidel_RB_2D, Poisson_Jacobi_2D, &
      Poisson_assemble_matrix_2D, Poisson_residual_2D
   implicit none

   private

   !> Number of dimensions and number of derivatives
   integer, parameter :: ndims = 2
   integer, parameter :: num_derivatives = 4
   integer, dimension(ndims*num_derivatives), parameter :: derivatives = [1,0,0,1,2,0,0,2] ! dx, dy, dxx, dyy

   ! Grid parameters
   integer, dimension(ndims), parameter :: grid_size = [64,64], processor_dims = [1,1]
   logical, dimension(ndims), parameter :: periods = [.false.,.false.]
   logical, parameter :: reorder = .true.
   real, dimension(ndims), parameter :: domain_begin = [0,0], domain_end = [20.0,10.0]
   integer, dimension(ndims), parameter :: stencil_sizes = 3
   integer, dimension(ndims), parameter :: uv_ghost_begin = [0,0], uv_ghost_end = [0,0]
   integer, dimension(ndims), parameter :: p_ghost_begin = [0,0], p_ghost_end = [0,0]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   ! Solver parameters
   integer, parameter :: t_steps = 5000
   real, parameter :: t0 = 0.0, t1 = 1.0, dt = (t1-t0)/t_steps

   ! Physical parameters
   real, parameter :: rho = 1.0, nu = 1.0, u_inflow = 100.0, cylinder_radius = 0.5, tol_cylinder = 1.0e-6
   real, dimension(ndims), parameter :: F = [0.0, 0.0], cylinder_center = [5.0,5.0]

   public :: flow_around_cylinder_main

contains

   !> Solve the 2D Navier-Stokes equation
   subroutine flow_around_cylinder_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      type(SolverParamsType):: solver
      type(comm_type) :: comm
      type(FDstencil_type) :: FDstencil
      type(block_type) :: u1_block, v1_block, p_block, &
         u2_block, v2_block, &
         u3_block, v3_block, &
         u4_block, v4_block
      type(ResultType) :: result

      integer :: ii, iounit
      real, dimension(4) :: result_array_with_timings

      write(0,*) "dt: ", dt

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, u1_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, v1_block)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, u2_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, v2_block)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, u3_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, v3_block)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, u4_block)
      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         uv_ghost_begin, uv_ghost_end, stencil_begin, stencil_end, 1, v4_block)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm, &
         p_ghost_begin, p_ghost_end, stencil_begin, stencil_end, 0, p_block)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil)

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      !> Assemble the matrix A
      call Poisson_assemble_matrix_2D(2, 3, p_block, FDstencil)

      ! Apply the matrix boundary conditions before decomposing the matrix
      call write_matrix_pressure_bc_2D(p_block, FDstencil)

      !> Decompose the matrix A into LU
      call LU_decomposition(p_block%direct_solver_matrix_ptr_2D, p_block%ipiv)

      ! Run the solver
      do ii = 1, t_steps/10
         call RK4(dt, comm, p_block, &
            u1_block, v1_block, &
            u2_block, v2_block, &
            u3_block, v3_block, &
            u4_block, v4_block, &
            FDstencil, solver, result)

         if(rank == MASTER_RANK) then
            print*, "Iteration: ", ii, "/", t_steps
            call print_resultType(result)

            if(result%converged == -1) then
               exit
            end if
         end if
      end do

      result_array_with_timings(2) = MPI_WTIME()

      result_array_with_timings(3) = result_array_with_timings(2) - result_array_with_timings(1)

      call all_reduce_mpi_wrapper(result_array_with_timings(3), result_array_with_timings(4), 1, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm%comm)

      ! Write out the result and timings from the master rank
      if(rank == MASTER_RANK) then
         write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(4)/world_size, " seconds"
      end if

      ! Write out our setup to a file. Just for debugging purposes
      call open_txt_file("output/output_from_", rank, iounit)

      !call print_cart_comm_type(comm_params, iounit)

      call print_block_type(u1_block, iounit)
      call print_block_type(v1_block, iounit)
      call print_block_type(p_block, iounit)

      !call print_finite_difference_stencil(FDstencil_params, iounit)

      call close_txt_file(iounit)

      ! Write out system solution to a file
      call write_block_data_to_file(u1_block%data_layout, "output/u_solution.dat", comm%comm, u1_block%matrix)
      call write_block_data_to_file(v1_block%data_layout, "output/v_solution.dat", comm%comm, v1_block%matrix)
      call write_block_data_to_file(p_block%data_layout, "output/p_solution.dat", comm%comm, p_block%matrix)

      ! Deallocate data

      call deallocate_cart_comm_type(comm)

      call deallocate_block_type(p_block)

      call deallocate_block_type(u1_block)
      call deallocate_block_type(v1_block)

      call deallocate_block_type(u2_block)
      call deallocate_block_type(v2_block)

      call deallocate_block_type(u3_block)
      call deallocate_block_type(v3_block)

      call deallocate_block_type(u4_block)
      call deallocate_block_type(v4_block)

      call deallocate_finite_difference_stencil(FDstencil)

   end subroutine flow_around_cylinder_main

   !> Subroutine to timestep using the RK4-method
   subroutine RK4(dt, comm, p_block, &
      u1_block, v1_block, &
      u2_block, v2_block, &
      u3_block, v3_block, &
      u4_block, v4_block, &
      FDstencil, solver, result)
      real, intent(in) :: dt
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: p_block, &
         u1_block, v1_block, &
         u2_block, v2_block, &
         u3_block, v3_block, &
         u4_block, v4_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver
      type(ResultType), intent(out) :: result

      !> RK4 time coefficients
      real, parameter :: RK1_time = 0, RK2_time = 1.0/2.0, RK3_time = 1.0/2.0, RK4_time = 1.0

      !> RK4 combination coefficients
      real, parameter :: RK1_coef = 1.0/6.0, RK2_coef = 1.0/3.0, RK3_coef = 1.0/3.0, RK4_coef = 1.0/6.0

      integer :: ii

      ! Use the RK4 method to find the intermediate velocity fields only!
      ! f_matrix_ptr is the rhs found from the intermediate velocity fields

      !> RK4 step 1
      call calculate_intermediate_velocity_2D(dt, u1_block, v1_block, FDstencil)
      ! uk1 = dt * rhs(u1) Done in calculate_intermediate_velocity_2D
      ! vk1 = dt * rhs(v1) Done in calculate_intermediate_velocity_2D
      call daxpy_to_vector(RK2_time, u1_block%matrix_ptr, u1_block%f_matrix_ptr, u2_block%matrix_ptr) ! u2 = u1 + 1/2 * uk1
      call daxpy_to_vector(RK2_time, v1_block%matrix_ptr, v1_block%f_matrix_ptr, v2_block%matrix_ptr) ! v2 = v1 + 1/2 * vk1
      call project_velocity(RK2_time, comm, u2_block, v2_block, p_block, FDstencil, solver, result) ! Project the velocity fields

      !> RK4 step 2
      call calculate_intermediate_velocity_2D(dt,u2_block, v2_block, FDstencil)
      ! uk2 = dt * rhs(u2) Done in calculate_intermediate_velocity_2D
      ! vk2 = dt * rhs(v2) Done in calculate_intermediate_velocity_2D
      call daxpy_to_vector(RK3_time, u1_block%matrix_ptr, u2_block%f_matrix_ptr, u3_block%matrix_ptr) ! u3 = u1 + 1/2 * uk2
      call daxpy_to_vector(RK3_time, v1_block%matrix_ptr, v2_block%f_matrix_ptr, v3_block%matrix_ptr) ! v3 = v1 + 1/2 * vk2
      call project_velocity(RK3_time, comm, u3_block, v3_block, p_block, FDstencil, solver, result) ! Project the velocity fields

      !> RK4 step 3
      call calculate_intermediate_velocity_2D(dt, u3_block, v3_block, FDstencil)
      ! uk3 = dt * rhs(u3) Done in calculate_intermediate_velocity_2D
      ! vk3 = dt * rhs(v3) Done in calculate_intermediate_velocity_2D
      call daxpy_to_vector(RK4_time, u1_block%matrix_ptr, u3_block%f_matrix_ptr, u4_block%matrix_ptr) ! u4 = u1 + 1 * uk3
      call daxpy_to_vector(RK4_time, v1_block%matrix_ptr, v3_block%f_matrix_ptr, v4_block%matrix_ptr) ! v4 = v1 + 1 * vk3
      call project_velocity(RK4_time, comm, u4_block, v4_block, p_block, FDstencil, solver, result) ! Project the velocity fields

      !> RK4 step 4
      call calculate_intermediate_velocity_2D(dt, u4_block, v4_block, FDstencil)
      ! uk4 = dt * rhs(u4) Done in calculate_intermediate_velocity_2D
      ! vk4 = dt * rhs(v4) Done in calculate_intermediate_velocity_2D

      !> Combine the results
      !$omp parallel do default(none) &
      !$omp shared(u1_block, u2_block, u3_block, u4_block, v1_block, v2_block, v3_block, v4_block) &
      !$omp private(ii)
      do ii = 1, size(u1_block%matrix_ptr)
         ! u_intermediate = u1 + 1/6 * (uk1 + 2*uk2 + 2*uk3 + uk4)
         u1_block%matrix_ptr(ii) = u1_block%matrix_ptr(ii) + (RK1_coef * u1_block%f_matrix_ptr(ii) + &
            RK2_coef * u2_block%f_matrix_ptr(ii) + RK3_coef * u3_block%f_matrix_ptr(ii) + RK4_coef * u4_block%f_matrix_ptr(ii))

         ! v_intermediate = v1 + 1/6 * (vk1 + 2*vk2 + 2*vk3 + vk4)
         v1_block%matrix_ptr(ii) = v1_block%matrix_ptr(ii) + (RK1_coef * v1_block%f_matrix_ptr(ii) + &
            RK2_coef * v2_block%f_matrix_ptr(ii) + RK3_coef * v3_block%f_matrix_ptr(ii) + RK4_coef * v4_block%f_matrix_ptr(ii))
      end do
      !$omp end parallel do

      ! Project the final velocity fields
      call project_velocity(dt, comm, u1_block, v1_block, p_block, FDstencil, solver, result)

   end subroutine RK4

   ! Solve the Navier-Stokes in 2D using the fractional step method.
   subroutine project_velocity(dt, comm, u_block, v_block, p_block, FDstencil, solver, result)
      real, intent(in) :: dt
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver
      type(ResultType), intent(out) :: result

      !> Write out the boundary conditions on the velocity fields
      call write_velocity_BC_2D(u_block%extended_block_begin_c, u_block%extended_block_end_c, &
         u_block%matrix_ptr_2D, v_block%matrix_ptr_2D)

      ! Calculate velocity divergence
      call calculate_velocity_divergence_2D(dt, u_block, v_block, p_block, FDstencil)

      ! Calculate pressure correction directly
      call set_bc_zero_2D(p_block%f_matrix_ptr_2D)
      !call write_rhs_pressure_bc_2D(p_block, FDstencil)
      call solve_LU_system(p_block%direct_solver_matrix_ptr_2D, p_block%f_matrix_ptr, p_block%ipiv)
      call copy_vector(p_block%f_matrix_ptr, p_block%matrix_ptr)

      ! Update the velocity fields and correct the pressure
      call update_velocity_fields_2D(dt, u_block, v_block, p_block, FDstencil)

      ! Write out the boundary conditions on the velocity fields
      call write_velocity_BC_2D(u_block%extended_block_begin_c, u_block%extended_block_end_c, &
         u_block%matrix_ptr_2D, v_block%matrix_ptr_2D)

   end subroutine project_velocity

   ! Write the boundary conditions on the velocity fields.
   subroutine write_velocity_BC_2D(extended_begin_c, extended_end_c, u_matrix, v_matrix)
      integer, dimension(:), intent(in) :: extended_begin_c, extended_end_c
      real, dimension(:,:), intent(inout) :: u_matrix, v_matrix

      integer :: ii, jj, on_boundary
      integer, dimension(2) :: local_indices
      real, dimension(2) :: point, normal

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(extended_begin_c, extended_end_c, u_matrix, v_matrix) &
      !$omp private(ii, jj, local_indices, point, normal, on_boundary)
      do jj = extended_begin_c(2)+1, extended_end_c(2)
         do ii = extended_begin_c(1)+1, extended_end_c(1)
            local_indices = [ii,jj]

            call calculate_point(2, [0,0], local_indices, domain_begin, (domain_end-domain_begin)/(grid_size-1), point)

            ! Left wall: u = inflow, v = 0
            if(local_indices(1) == 1) then
               u_matrix(local_indices(1),local_indices(2)) = u_inflow*(1.0-((point(2)-domain_end(2)/2)/domain_end(2)/2)**2)
               v_matrix(local_indices(1),local_indices(2)) = 0.0
            end if

            ! Right wall: u_x = 0, v_x = 0
            if(local_indices(1) == extended_end_c(1)) then
               u_matrix(local_indices(1),local_indices(2)) = u_matrix(local_indices(1)-1,local_indices(2))
               v_matrix(local_indices(1),local_indices(2)) = v_matrix(local_indices(1)-1,local_indices(2))
            end if

            ! Top wall: u = 0, v = 0
            if (local_indices(2) == 1) then
               u_matrix(local_indices(1),local_indices(2)) = u_matrix(local_indices(1),local_indices(2)+1)
               v_matrix(local_indices(1),local_indices(2)) = 0.0
            end if

            ! Bottom wall: u = 0, v = 0
            if(local_indices(2) == extended_end_c(2)) then
               u_matrix(local_indices(1),local_indices(2)) = u_matrix(local_indices(1),local_indices(2)-1)
               v_matrix(local_indices(1),local_indices(2)) = 0.0
            end if

            ! Cylinder wall and inside the cylinder: u = 0, v = 0
            call on_cylinder_boundary(point, cylinder_center, cylinder_radius, cylinder_radius, on_boundary, normal)
            if(on_boundary == 1) then
               u_matrix(local_indices(1),local_indices(2)) = 0.0
               v_matrix(local_indices(1),local_indices(2)) = 0.0
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_velocity_BC_2D

   ! Calculate rhs for the 2D Navier-Stokes equation
   subroutine calculate_intermediate_velocity_2D(dt, u_block, v_block, FDstencil)
      real, intent(in) :: dt
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
      !$omp shared(dt, u_block, v_block, FDstencil) &
      !$omp private(ii, jj, uv_local_indices, u, v, u_rhs, v_rhs, coefficients, alpha, beta, dx, dy, dxx, dyy, combined_stencil, &
      !$omp dx_2D, dy_2D, dxx_2D, dyy_2D, combined_stencil_2D)
      do jj = u_block%block_begin_c(2)+1, u_block%block_end_c(2)
         do ii = u_block%block_begin_c(1)+1, u_block%block_end_c(1)
            uv_local_indices = [ii,jj]

            u = u_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2))
            v = v_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2))

            call get_coefficients_wrapper(FDstencil, [1,1], u_block%extended_block_dims, &
               uv_local_indices, alpha, beta, coefficients)

            dx => coefficients(1:FDstencil%num_stencil_elements)
            dy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)
            dxx => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
            dyy => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

            combined_stencil = -(u*dx + v*dy) + nu * (dxx + dyy) + F(1)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencil, combined_stencil_2D)
            call apply_FDstencil_2D(combined_stencil_2D, u_block%matrix_ptr_2D, uv_local_indices, alpha, beta, u_rhs)

            combined_stencil = -(u*dx + v*dy) + nu * (dxx + dyy) + F(2)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencil, combined_stencil_2D)
            call apply_FDstencil_2D(combined_stencil_2D, v_block%matrix_ptr_2D, uv_local_indices, alpha, beta, v_rhs)

            u_block%f_matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2)) = dt * u_rhs ! uk = dt * rhs(u)
            v_block%f_matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2)) = dt * v_rhs ! vk = dt * rhs(v)

         end do
      end do

      !$omp end parallel do

   end subroutine calculate_intermediate_velocity_2D

   !> Calculate velocity divergence
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

            dx => coefficients(1:FDstencil%num_stencil_elements)
            dy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dy, dy_2D)

            call apply_FDstencil_2D(dx_2D, u_block%matrix_ptr_2D, uv_local_indices, alpha, beta, u_x)
            call apply_FDstencil_2D(dy_2D, v_block%matrix_ptr_2D, uv_local_indices, alpha, beta, v_y)

            p_block%f_matrix_ptr_2D(p_local_indices(1), p_local_indices(2)) = (rho/dt) * (u_x + v_y)

         end do
      end do

      !$omp end parallel do

   end subroutine calculate_velocity_divergence_2D

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

            dx => coefficients(1:FDstencil%num_stencil_elements)
            dy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dy, dy_2D)

            call apply_FDstencil_2D(dx_2D, p_block%matrix_ptr_2D, p_local_indices, alpha, beta, p_x)
            call apply_FDstencil_2D(dy_2D, p_block%matrix_ptr_2D, p_local_indices, alpha, beta, p_y)

            u_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2)) = &
               u_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2))  - (dt/rho) * p_x
            v_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2))  = &
               v_block%matrix_ptr_2D(uv_local_indices(1),uv_local_indices(2))  - (dt/rho) * p_y

         end do
      end do

      !$omp end parallel do

   end subroutine update_velocity_fields_2D

   !> Write the Dirichlet and Neumann boundary conditions to the matrix
   subroutine write_matrix_pressure_bc_2D(p_block, FDstencil)
      type(block_type), intent(inout) :: p_block
      type(FDstencil_type), intent(inout) :: FDstencil

      integer :: ii, jj, diag, on_boundary
      integer, dimension(ndims) :: local_indices, global_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients
      real, contiguous, dimension(:), pointer :: dx, dy
      real, dimension(2) :: point, normal

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(p_block, FDstencil) &
      !$omp private(ii, jj, local_indices, global_indices, diag, alpha, beta, coefficients, dx, dy, point, on_boundary, normal)
      do jj = p_block%extended_block_begin_c(2)+1, p_block%extended_block_end_c(2)
         do ii = p_block%extended_block_begin_c(1)+1, p_block%extended_block_end_c(1)
            local_indices = [ii,jj]
            global_indices = p_block%extended_global_begin_c + local_indices

            call calculate_point(2, [0,0], local_indices, domain_begin, (domain_end-domain_begin)/(grid_size-1), point)

            call IDX_XD(p_block%ndims, p_block%extended_block_dims, local_indices, diag)

            call get_coefficients_wrapper(FDstencil, [1,1], p_block%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dx => coefficients(1:FDstencil%num_stencil_elements)
            dy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            ! Left wall, p_x = 0 (Neumann condition)
            if(local_indices(1) == 1) then
               call set_matrix_coefficients(p_block%ndims, FDstencil%stencil_sizes, p_block%extended_block_dims, &
                  dx, p_block%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)
            end if

            ! Right wall, p = 0 Dirichlet condition
            if(local_indices(1) == p_block%extended_block_end_c(1)) then
               p_block%direct_solver_matrix_ptr_2D(diag,:) = 0.0
               p_block%direct_solver_matrix_ptr_2D(diag,diag) = 1.0
            end if

            ! Top wall, p_y = 0 (Neumann condition)
            if(local_indices(2) == 1) then
               call set_matrix_coefficients(p_block%ndims, FDstencil%stencil_sizes, p_block%extended_block_dims, &
                  dy, p_block%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)
            end if

            ! Bottom wall, p_y = 0 (Neumann condition)
            if(local_indices(2) == p_block%extended_block_end_c(2)) then
               call set_matrix_coefficients(p_block%ndims, FDstencil%stencil_sizes, p_block%extended_block_dims, &
                  dy, p_block%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)
            end if

            ! Set the pressure to zero on the whole cylinder
            call on_cylinder_boundary(point, cylinder_center, cylinder_radius, cylinder_radius, on_boundary, normal)
            if(on_boundary == 1) then
               !p_block%direct_solver_matrix_ptr_2D(diag,:) = 0.0
               !p_block%direct_solver_matrix_ptr_2D(diag,diag) = 1.0
            end if

            ! Set the normal*grad(p) = 0 on the cylinder boundary
            call on_cylinder_boundary(point, cylinder_center, cylinder_radius, tol_cylinder, on_boundary, normal)
            if(on_boundary == 1) then
               call set_matrix_coefficients(p_block%ndims, FDstencil%stencil_sizes, p_block%extended_block_dims, &
                  normal(1)*dx+normal(2)*dy, p_block%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_matrix_pressure_bc_2D

   !> Write the rhs boundary conditions to the vector to be solve by the matrix
   subroutine write_rhs_pressure_bc_2D(p_block, FD_stencil)
      type(block_type), intent(inout) :: p_block
      type(FDstencil_type), intent(inout) :: FD_stencil

      integer :: ii, jj, on_boundary
      integer, dimension(ndims) :: local_indices, alpha, beta
      real, dimension(2) :: point, normal

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(p_block, FD_stencil) &
      !$omp private(ii, jj, local_indices, point, on_boundary, normal)
      do jj = p_block%extended_block_begin_c(2)+1, p_block%extended_block_end_c(2)
         do ii = p_block%extended_block_begin_c(1)+1, p_block%extended_block_end_c(1)
            local_indices = [ii,jj]

            call calculate_point(2, [0,0], local_indices, domain_begin, (domain_end-domain_begin)/(grid_size-1), point)

            call on_cylinder_boundary(point, cylinder_center, cylinder_radius, cylinder_radius, on_boundary, normal)
            if(on_boundary == 1) then
               p_block%f_matrix_ptr_2D(local_indices(1),local_indices(2)) = 0.0
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_rhs_pressure_bc_2D

   ! Check if a point is on the boundary of the cylinder
   pure subroutine on_cylinder_boundary(point, center, R, tol, on_boundary, normal)
      real, dimension(:), intent(in) :: point, center
      real, intent(in) :: R, tol
      integer, intent(out) :: on_boundary
      real, dimension(:), intent(out) :: normal

      real :: distance

      ! Find the distance from the cylinder center using built in vector operations
      distance = sqrt(sum((point - center)**2))

      if(abs(distance - R) <= tol) then
         on_boundary = 1
         normal = (point - center) / distance
      else
         on_boundary = 0
         normal = 0.0
      end if

   end subroutine on_cylinder_boundary

end module flow_around_cyliner_module

