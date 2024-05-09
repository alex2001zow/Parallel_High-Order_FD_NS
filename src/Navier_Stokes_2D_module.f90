module Navier_Stokes_2D_module
   use constants_module, only: pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type, print_cartesian_grid
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, &
      apply_FDstencil, update_value_from_stencil, calculate_scaled_coefficients, get_FD_coefficients_from_index
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type
   use functions_module, only: FunctionPtrType, FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers
   use initialization_module, only: write_function_to_block, write_initial_condition_and_boundary

   use solver_module, only: SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, check_convergence
   implicit none

   private

   integer, parameter :: N_iterations = 200
   integer, parameter :: N_Pressure_Poisson_iterations = 5000
   real, parameter :: Pressure_Poisson_tol = 1e-6
   real, parameter :: system_rho = 1.0, system_nu = 0.1, system_u_lid = 10.0, system_dt = 0.0001
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

      grid_size = 128
      processor_dims = 1
      domain_begin = 0
      domain_end = 1
      stencil_sizes = 3

      ! Set the solver parameters tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(Pressure_Poisson_tol, N_Pressure_Poisson_iterations, 2, solver_params)

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

      integer :: iounit, ii, jj, global_index, converged, iter
      real, dimension(4) :: result_array_with_timings

      real :: dt

      type(comm_type) :: comm_params
      type(FDstencil_type) :: FDstencil_params
      type(block_type) :: u_block_params, v_block_params, p_block_params, test_block_params
      type(FunctionPtrType) :: u_test_2D_func

      u_test_2D_func%output_size = 1
      u_test_2D_func%func => u_test_2D

      call create_cart_comm_type(ndims, processor_dims, rank, world_size, comm_params)

      call create_block_type(ndims, 1, 1, stencil_sizes * 0, stencil_sizes, domain_begin, domain_end, grid_size, &
         comm_params, u_block_params)

      call create_block_type(ndims, 1, 1, stencil_sizes * 0, stencil_sizes, domain_begin, domain_end, grid_size, &
         comm_params, v_block_params)

      call create_block_type(ndims, 1, 1, stencil_sizes * 0, stencil_sizes, &
         domain_begin, domain_end, grid_size, comm_params, p_block_params) ! Extra points for pressure due to neumann BC

      !call sleeper_function(1)

      ! call create_block_type(ndims, 1, 1, grid_size * 0, domain_begin, domain_end, grid_size, comm_params, test_block_params)

      ! call write_function_to_block(test_block_params%ndims, 1,test_block_params%domain_begin, test_block_params%domain_end, &
      !    test_block_params%size, test_block_params%begin, test_block_params%size, test_block_params%matrix, &
      !    test_block_params%dx, u_test_2D_func)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      !! Scale the coefficients by dx
      call calculate_scaled_coefficients(ndims, u_block_params%dx, FDstencil_params)

      ! Initialize the block
      u_block_params%matrix = 0.0
      v_block_params%matrix = 0.0
      p_block_params%matrix = 0.0

      !call write_velocity_BC(ndims, u_block_params%extended_local_size, u_block_params%matrix, v_block_params%matrix)
      !call write_pressure_BC(ndims, FDstencil_params, p_block_params%extended_local_size, p_block_params%matrix)

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      dt = system_dt
      ! Run the solver
      do ii = 1, N_iterations
         call NS_2D_one_timestep(comm_params, u_block_params, v_block_params, p_block_params, &
            FDstencil_params, solver_params, system_rho, system_nu, dt, converged, iter)
         write(*,*) "Iteration: ", ii, "/", N_iterations, " Converged: ", converged, " Iter: ", iter

         if(converged == -1) then
            exit
         end if

      end do

      !call test_fd_method(test_block_params, FDstencil_params)

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

      call print_finite_difference_stencil(ndims, FDstencil_params, iounit)

      call print_block_type(ndims, u_block_params, iounit)
      call print_block_type(ndims, v_block_params, iounit)
      call print_block_type(ndims, p_block_params, iounit)

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

   end subroutine Navier_Stokes_2D_analytical_wrapper

   ! Solve the Navier-Stokes in 2D using the fractional step method.
   subroutine NS_2D_one_timestep(comm, u_block, v_block, p_block, FDstencil, solver_params_in, rho, nu, dt, converged, it)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: u_block, v_block, p_block
      type(FDstencil_type), target, intent(inout) :: FDstencil
      type(SolverParamsType), intent(in) :: solver_params_in
      real, intent(in) :: rho, nu, dt
      integer, intent(inout) :: converged
      integer, intent(inout) :: it

      real, dimension(product(u_block%extended_local_size)) :: u_rhs, v_rhs, u_star, v_star, u_x_v_y, p_temp

      integer :: ndims, num_elements

      real, dimension(4) :: norm_array

      ndims = comm%ndims
      num_elements = 1
      norm_array = [1e12,1e18,1e36,1e48]

      u_rhs = 0
      v_rhs = 0
      u_star = 0
      v_star = 0
      u_x_v_y = 0

      ! Calculate intermediate velocity fields
      call compute_rhs(u_block%ndims, u_block%extended_local_size, u_block%extended_block_begin_c+1, num_elements, &
         u_block%matrix, v_block%matrix, FDstencil, nu, u_rhs, v_rhs)
      u_star = u_block%matrix + dt * u_rhs
      v_star = v_block%matrix + dt * v_rhs

      !call compute_rhs(u_block%ndims, u_block%extended_block_dims, u_block%extended_block_begin_c+1, num_elements, &
      !   u_star, v_star, FDstencil, u_rhs, v_rhs)
      !u_star = u_block%matrix + 0.5*dt*(u_rhs + (u_star - u_block%matrix)/dt)
      !v_star = v_block%matrix + 0.5*dt*(v_rhs + (v_star - v_block%matrix)/dt)

      call write_velocity_BC(u_block%ndims, u_block%extended_block_dims, u_star, v_star)

      ! Calculate velocity divergence
      call calculate_velocity_divergence(u_block%ndims, u_block%extended_block_dims, u_block%block_begin_c+1, num_elements, &
         u_star, v_star, FDstencil, u_x_v_y)

      p_temp = p_block%matrix

      ! Calculate pressure correction
      call solve_poisson_problem(comm%comm, p_block%ndims, p_block%extended_block_dims, u_block%extended_block_dims, &
         p_block%extended_block_begin_c+1, num_elements, &
         p_block%matrix, p_temp, u_x_v_y, FDstencil, solver_params_in%tol, solver_params_in%max_iter, &
         norm_array, rho, dt, converged, it)

      ! Update the velocity fields and correct the pressure
      call update_velocity_fields(ndims, u_block, v_block, p_block, FDstencil, p_block%matrix, u_star, v_star, rho, dt)

      ! Write out the boundary conditions on the velocity fields
      call write_velocity_BC(u_block%ndims, u_block%extended_block_dims, u_star, v_star)

      !call normalize_pressure_matrix(p_block%ndims, p_block%size, p_block%matrix)

   end subroutine NS_2D_one_timestep

   ! Write the lid cavity boundary conditions on the velocity fields.
   subroutine write_velocity_BC(ndims, dims, u_matrix, v_matrix)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: dims
      real, dimension(:), intent(inout) :: u_matrix, v_matrix

      integer :: ii, jj, global_index

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(ndims, dims, u_matrix, v_matrix) &
      !$omp private(ii, jj, global_index)
      do ii = 1, dims(1)
         do jj = 1, dims(2)
            call IDX_XD(ndims, dims, [ii,jj], global_index)

            if(ii == 1 .or. ii == dims(1) .or. jj == 1 .or. jj == dims(2)) then
               u_matrix(global_index) = 0.0
               v_matrix(global_index) = 0.0
            end if

            if(ii == 1) then
               u_matrix(global_index) = system_u_lid
            end if
         end do
      end do

      !$omp end parallel do

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
      integer, dimension(ndims) :: p_index, alpha
      real, dimension(:), pointer :: coefficients, dfx, dfy
      real, dimension(product(FDstencil%stencil_sizes)) :: combined_stencil
      real :: old_val, new_val

      combined_stencil = 0

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(ndims, p_dims, FDstencil, p_matrix) &
      !$omp private(ii, jj, p_index, global_index, alpha, coefficients, dfx, dfy, old_val, new_val)
      do ii = 1, p_dims(1)
         do jj = 1, p_dims(2)
            p_index = [ii, jj]
            call IDX_XD(ndims, p_dims, p_index, global_index)

            call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               [1, 1], p_dims, p_index, FDstencil%scaled_stencil_coefficients, alpha, coefficients)

            dfx => coefficients(1:FDstencil%num_stencil_elements)
            dfy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            old_val = 0.0!p_matrix(global_index)

            ! Left and right walls: u_x = 0 (Neumann condition)
            ! if (jj == 1 .or. jj == p_dims(2)) then
            !    call update_value_from_stencil(ndims, 1, 0, FDstencil%stencil_sizes, alpha, dfx, &
            !       p_dims, p_index, p_matrix, old_val, new_val)

            !    p_matrix(global_index) = new_val
            ! end if

            ! Bottom wall: u_y = 0 (Neumann condition)
            ! if (ii == p_dims(1)) then
            !    call update_value_from_stencil(ndims, 1, 0, FDstencil%stencil_sizes, alpha, dfy, &
            !       p_dims, p_index, p_matrix, old_val, new_val)

            !    p_matrix(global_index) = new_val
            ! end if

            ! Corners: u_x + u_y = 0 (Neumann condition)
            ! if ((ii == 1 .or. ii == p_dims(1)) .and. (jj == 1 .or. jj == p_dims(2))) then
            !    combined_stencil = (dfx+dfy)
            !    call update_value_from_stencil(ndims, 1, 0, FDstencil%stencil_sizes, alpha, combined_stencil, &
            !       p_dims, p_index, p_matrix, old_val, new_val)

            !    p_matrix(global_index) = new_val
            ! end if

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
      integer, dimension(ndims) :: index, alphas
      integer, dimension(:), pointer :: stencil_size
      real, dimension(:), pointer :: coefficients, dfx, dfy, dfxx, dfyy
      real :: u, v, u_x, u_y, v_x, v_y, u_xx, u_yy, v_xx, v_yy, f_u, f_v

      stencil_size => FDstencil%stencil_sizes

      ! Perhaps combine the stencils into a single one. I think it is doable.

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(dims, u_matrix, v_matrix, FDstencil, stencil_size, nu, u_rhs, v_rhs, num_elements, ndims, start_dims) &
      !$omp private(ii, jj, index, global_index, alphas, coefficients, dfx, dfy, dfxx, dfyy, &
      !$omp u, v, u_x, u_y, u_xx, u_yy, v_x, v_y, v_xx, v_yy, f_u, f_v)
      do ii = 2, dims(1)-1
         do jj = 2, dims(2)-1
            index = [ii,jj]
            call IDX_XD(ndims, dims, index, global_index)

            call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               start_dims, dims, index, FDstencil%scaled_stencil_coefficients, alphas, coefficients)

            dfx => coefficients(1:FDstencil%num_stencil_elements)
            dfy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)
            dfxx => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
            dfyy => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfx, dims, index, u_matrix, u_x)
            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfxx, dims, index, u_matrix, u_xx)
            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfy, dims, index, u_matrix, u_y)
            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfyy, dims, index, u_matrix, u_yy)

            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfx, dims, index, v_matrix, v_x)
            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfxx, dims, index, v_matrix, v_xx)
            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfy, dims, index, v_matrix, v_y)
            call apply_FDstencil(ndims, num_elements, 0, stencil_size, alphas, dfyy, dims, index, v_matrix, v_yy)

            u = u_matrix(global_index)
            v = v_matrix(global_index)

            f_u = 0
            f_v = 0

            u_rhs(global_index) = -(u*u_x + v * u_y) + nu * (u_xx + u_yy) + f_u
            v_rhs(global_index) = -(u*v_x + v * v_y) + nu * (v_xx + v_yy) + f_v

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
      integer, dimension(ndims) :: index, alphas
      integer,dimension(:), pointer :: stencil_size
      real, dimension(:), pointer :: coefficients, dfx, dfy
      real :: u_x, v_y

      stencil_size => FDstencil%stencil_sizes

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(dims, ndims, num_elements, start_dims, FDstencil, stencil_size, u_matrix, v_matrix, u_x_v_y) &
      !$omp private(ii, jj, index, global_index, alphas, coefficients, dfx, dfy, u_x, v_y)
      do ii = 2, dims(1)-1
         do jj = 2, dims(2)-1
            index = [ii,jj]
            call IDX_XD(ndims, dims, index, global_index)

            call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               start_dims, dims, index, FDstencil%scaled_stencil_coefficients, alphas, coefficients)

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
      integer, dimension(ndims) :: p_index, alphas
      integer, dimension(:), pointer :: stencil_size

      real, dimension(:), pointer :: coefficients, dfxx, dfyy, p_ptr, p_temp_ptr
      real, dimension(product(FDstencil%stencil_sizes)) :: combined_stencil
      real :: local_norm, old_val, new_val, f_val

      stencil_size => FDstencil%stencil_sizes

      p_ptr => p_matrix
      p_temp_ptr => p_temp_matrix

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < max_iter)

         local_norm = 0.0

         !$omp parallel do collapse(2) default(none) &
         !$omp shared(p_dims, ndims, rho, dt, u_x_v_y, p_ptr, p_temp_ptr, FDstencil, stencil_size, num_elements, p_start_dims) &
         !$omp private(ii, jj, p_index, global_index, old_val, alphas, coefficients, dfxx, dfyy, combined_stencil, f_val, new_val) &
         !$omp reduction(+:local_norm)
         do ii = 2, p_dims(1)-1
            do jj = 2, p_dims(2)-1
               p_index = [ii,jj]
               call IDX_XD(ndims, p_dims, p_index, global_index)

               old_val = p_ptr(global_index)

               call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
                  p_start_dims, p_dims, p_index, FDstencil%scaled_stencil_coefficients, alphas, coefficients)

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

         call check_convergence(comm, tol, 1.0/product(p_dims), norm_array, converged)
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
      integer, dimension(ndims) :: p_start_dims, p_index, alphas
      integer, dimension(:), pointer :: stencil_size
      real, dimension(:), pointer :: coefficients, dfx, dfy
      real :: p_corr_x, p_corr_y

      stencil_size => FDstencil%stencil_sizes
      p_start_dims = p_block%block_begin_c+1

      ! Update the velocity fields

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(p_block, ndims, num_elements, p_start_dims, FDstencil, stencil_size, &
      !$omp p_corr_array, u_block, u_star, v_block, v_star, dt, rho) &
      !$omp private(ii, jj, p_index, p_global_index, alphas, coefficients, dfx, dfy, p_corr_x, p_corr_y)
      do ii = 1,p_block%extended_local_size(1)
         do jj = 1,p_block%extended_local_size(2)
            p_index = [ii,jj]
            call IDX_XD(ndims, p_block%extended_local_size, p_index, p_global_index)

            call get_FD_coefficients_from_index(ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               p_start_dims, p_block%extended_local_size, p_index, FDstencil%scaled_stencil_coefficients, alphas, coefficients)

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

   pure subroutine normalize_pressure_matrix(ndims, dims, matrix)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: dims
      real, dimension(:), intent(inout) :: matrix

      ! Normalize the pressure field
      matrix = matrix - (sum(matrix) / product(dims))

   end subroutine normalize_pressure_matrix

   subroutine test_fd_method(test_block, FDstencil)
      type(block_type), intent(inout) :: test_block
      type(FDstencil_type), target, intent(inout) :: FDstencil

      integer :: ndims, ii, jj, global_index, num_elements

      real, dimension(:), pointer :: coefficients
      real, dimension(:), pointer :: dfx, dfy, dfxx, dfyy
      integer, dimension(test_block%ndims) :: stencil_size, alphas, start_dims, index

      real :: u_x, u_y, u_xx, u_yy

      real, dimension(1) :: u_x_test, u_y_test, u_xx_test, u_yy_test

      real, dimension(test_block%ndims) :: global_domain_begin, global_domain_end, dx
      integer, dimension(test_block%ndims) :: global_domain_size, global_begin_indices

      real :: error

      ndims = test_block%ndims

      global_domain_begin = test_block%domain_begin
      global_domain_end = test_block%domain_end

      global_domain_size = test_block%extended_local_size
      global_begin_indices = 1

      dx = test_block%dx

      num_elements = 1

      start_dims = 1

      stencil_size = FDstencil%stencil_sizes

      error = 0.0
      do ii = 1, test_block%extended_local_size(1)
         do jj = 1, test_block%extended_local_size(2)
            index = [ii,jj]
            call IDX_XD(test_block%ndims, test_block%extended_local_size, index, global_index)

            call get_FD_coefficients_from_index(test_block%ndims, FDstencil%num_derivatives, FDstencil%stencil_sizes, &
               start_dims, test_block%extended_local_size, index, FDstencil%scaled_stencil_coefficients, alphas, coefficients)

            dfx => coefficients(1:FDstencil%num_stencil_elements)
            dfy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)
            dfxx => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
            dfyy => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)

            call apply_FDstencil(test_block%ndims, num_elements, 0, stencil_size, alphas, dfx, &
               test_block%extended_local_size, index, test_block%matrix, u_x)
            call apply_FDstencil(test_block%ndims, num_elements, 0, stencil_size, alphas, dfy, &
               test_block%extended_local_size, index, test_block%matrix, u_y)
            call apply_FDstencil(test_block%ndims, num_elements, 0, stencil_size, alphas, dfxx, &
               test_block%extended_local_size, index, test_block%matrix, u_xx)
            call apply_FDstencil(test_block%ndims, num_elements, 0, stencil_size, alphas, dfyy, &
               test_block%extended_local_size, index, test_block%matrix, u_yy)

            call ux_test_2D(ndims, num_elements, global_begin_indices, index, &
               global_domain_begin, global_domain_end, global_domain_size, dx, u_x_test)
            call uy_test_2D(ndims, num_elements, global_begin_indices, index, &
               global_domain_begin, global_domain_end, global_domain_size, dx, u_y_test)
            call uxx_test_2D(ndims, num_elements, global_begin_indices, index, &
               global_domain_begin, global_domain_end, global_domain_size, dx, u_xx_test)
            call uyy_test_2D(ndims, num_elements, global_begin_indices, index, &
               global_domain_begin, global_domain_end, global_domain_size, dx, u_yy_test)

            error = error + abs(u_x - u_x_test(1)) + abs(u_y - u_y_test(1)) + abs(u_xx - u_xx_test(1)) + abs(u_yy - u_yy_test(1))

            !write(*,*) " u_x Diff: ", abs(u_x - u_x_test(1)), " u_y Diff: ", abs(u_y - u_y_test(1)), &
            !  " u_xx Diff: ", abs(u_xx - u_xx_test(1)), " u_yy Diff: ", abs(u_yy - u_yy_test(1))

         end do
      end do

      write(*,*) "Error: ", (error/4.0) * (1.0/product(test_block%extended_local_size))

   end subroutine test_fd_method

   pure subroutine u_test_2D(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      real, dimension(ndims) :: point
      real :: x, y

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      x = point(1)
      y = point(2)

      func_val(1) = 3*x*x + 2*x*y + y*y*y - 4*y + 7

   end subroutine u_test_2D

   pure subroutine ux_test_2D(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      real, dimension(ndims) :: point
      real :: x, y

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      x = point(1)
      y = point(2)

      func_val(1) = 6*x + 2*y

   end subroutine ux_test_2D

   pure subroutine uy_test_2D(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      real, dimension(ndims) :: point
      real :: x, y

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      x = point(1)
      y = point(2)

      func_val(1) = 2*x + 3*y*y - 4

   end subroutine uy_test_2D

   pure subroutine uxx_test_2D(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      real, dimension(ndims) :: point
      real :: x, y

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      x = point(1)
      y = point(2)

      func_val(1) = 6

   end subroutine uxx_test_2D

   pure subroutine uyy_test_2D(ndims, num_elements, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims, num_elements
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, dimension(num_elements), intent(inout) :: func_val

      real, dimension(ndims) :: point
      real :: x, y

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      x = point(1)
      y = point(2)

      func_val(1) = 6*y

   end subroutine uyy_test_2D

end module Navier_Stokes_2D_module
