module test_Poisson_module
   use constants_module, only : pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, apply_FDstencil, update_value_from_stencil, &
      calculate_scaled_coefficients, apply_FDstencil_2D, update_value_from_stencil_2D, &
      get_coefficients_wrapper, set_matrix_coefficients
   use block_module, only: block_type, create_block_type, deallocate_block_type, sendrecv_data_neighbors, print_block_type
   use functions_module, only: calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers, &
      swap_pointers_2D, reshape_real_1D_to_2D, print_real_2D_array, calculate_dt_from_CFL, &
      LU_decomposition, solve_LU_system, daxpy_wrapper, daxpy_to_vector

   use solver_module, only: SolverParamsType, set_SolverParamsType, &
      ResultType, set_ResultType, print_ResultType, check_convergence
   use multigrid_module, only: full_weighing_restriction_2D, bilinear_prolongation_2D, apply_correction
   use Poisson_functions_module, only: Poisson_Gauss_Seidel_2D, Poisson_Gauss_Seidel_RB_2D, &
      Poisson_Jacobi_2D, Poisson_assemble_matrix_2D, Poisson_residual_2D

   implicit none

   private

   !> Number of dimensions and number of derivatives
   integer, parameter :: ndims = 2, num_derivatives = 2
   integer, dimension(ndims*num_derivatives), parameter :: derivatives = [2,0,0,2] ! dxx, dyy

   !> Grid parameters
   integer, dimension(ndims), parameter :: grid_size = [64,64], processor_dims = [1,1]
   logical, dimension(ndims), parameter :: periods = [.false., .false.]
   logical, parameter :: reorder = .true.
   real, dimension(ndims), parameter :: domain_begin = [-10.0*pi,-10.0*pi], domain_end = [10.0*pi,10.0*pi]
   integer, dimension(ndims), parameter :: stencil_sizes = 5
   integer, dimension(ndims), parameter :: ghost_begin = [0,0], ghost_end = [0,0]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   ! Time parameters
   real, parameter :: t_0 = 0.0, t_1 = 1.0
   real :: current_time = t_0
   real :: t_steps, dt

   !> Solver parameters
   integer, parameter :: solve_iterativly = 0, Jacobi_or_GS = 1
   real, parameter :: tol = (1e-12)**2, div_tol = 1e-1, omega = 0.8
   integer, parameter :: max_iter = 10 * (2.0 - omega), multigrid_max_level = 4
   type(SolverParamsType) :: solver_params

   public :: Poisson_main

contains

   !> The main Poisson subroutine to solve the system
   subroutine Poisson_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      type(comm_type) :: comm_params
      type(FDstencil_type) :: FDstencil_params
      type(block_type) :: data_block
      type(ResultType) :: result

      integer :: iounit, ii
      real, dimension(4) :: result_array_with_timings

      !call sleeper_function(1)

      ! Set the solver parameters tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(tol, div_tol, max_iter, Jacobi_or_GS, solver_params)

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin, ghost_end, stencil_begin, stencil_end, solve_iterativly, data_block)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      !call calculate_dt_from_CFL(1.0, [1.0,1.0], data_block%extended_grid_dx, dt)
      !t_steps = (t_end - t_begin) / dt + 1
      t_steps = 1
      dt = (t_1 - t_0) / t_steps

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      if(solve_iterativly /= 1) then
         !> Assemble the matrix A
         call Poisson_assemble_matrix_2D(0, 1, data_block, FDstencil_params)

         ! Apply the Dirichlet boundary conditions to the matrix before decomposing it
         call write_matrix_dirichlet_bc_2D(data_block)

         !> Decompose the matrix A into LU
         call LU_decomposition(data_block%direct_solver_matrix_ptr_2D,data_block%ipiv)

         !> Write the initial condition to the system
         call write_solution_2D(data_block, t_0)
      end if

      ! Run the time loop
      do ii = 1, t_steps
         current_time = current_time + dt
         call one_timestep(data_block, FDstencil_params, comm_params, result)

         if(rank == MASTER_RANK) then
            print*, "Time step: ", ii, " Time: ", current_time
         end if

         call print_ResultType(result)

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
      !call open_txt_file("output/output_from_", rank, iounit)

      !call print_cart_comm_type(comm_params, iounit)

      !call print_finite_difference_stencil(FDstencil_params, iounit)

      !call print_block_type(data_block, iounit)

      !call close_txt_file(iounit)

      ! Write out system solution to a file
      if(solve_iterativly == 1) then
         call write_block_data_to_file(data_block%data_layout, "output/system_solution.dat", &
            comm_params%comm, data_block%matrix)
         call write_block_data_to_file(data_block%data_layout, "output/system_rhs.dat", &
            comm_params%comm, data_block%f_matrix)
         call write_block_data_to_file(data_block%data_layout, "output/system_residual.dat", &
            comm_params%comm, data_block%residual_matrix)
      else
         call write_block_data_to_file(data_block%data_layout, "output/system_solution.dat", &
            comm_params%comm, data_block%f_matrix)
         call write_block_data_to_file(data_block%data_layout, "output/system_rhs.dat", &
            comm_params%comm, data_block%matrix_ptr)
         call write_block_data_to_file(data_block%data_layout, "output/system_residual.dat", &
            comm_params%comm, data_block%f_matrix)
      end if

      ! Deallocate data

      call deallocate_cart_comm_type(comm_params)

      call deallocate_block_type(data_block)

      call deallocate_finite_difference_stencil(FDstencil_params)

   end subroutine Poisson_main

   subroutine one_timestep(data_block, FDstencil_params, comm_params, result)
      type(block_type), intent(inout) :: data_block
      type(FDstencil_type), intent(inout) :: FDstencil_params
      type(comm_type), intent(in) :: comm_params
      type(ResultType), intent(out) :: result

      !> Calculate the rhs of the system
      call write_rhs_2D(data_block, current_time)

      if(solve_iterativly == 1) then
         ! Apply the Dirichlet boundary conditions to the system solution
         call write_dirichlet_bc_2D(comm_params, data_block)

         ! Call the multigrid solver
         call w_cycle(solver_params, comm_params, data_block, FDstencil_params, result, multigrid_max_level, 1)

         ! Find the new solution at the next time step
         call daxpy_wrapper((dt/-2.0), data_block%matrix_ptr, data_block%new_matrix_ptr)

         ! Swap the pointers
         call swap_pointers(data_block%matrix_ptr, data_block%new_matrix_ptr)
      else
         ! Apply the Dirichlet boundary conditions to the rhs
         call write_rhs_dirichlet_bc_2D(data_block)

         ! Solve the linear system using the LU decomposition
         call solve_LU_system(data_block%direct_solver_matrix_ptr_2D, data_block%f_matrix_ptr, data_block%ipiv)

         ! Find the new solution at the next time step
         call daxpy_wrapper((dt/-2.0),data_block%f_matrix_ptr,data_block%new_matrix_ptr)

      end if

   end subroutine one_timestep

   recursive subroutine w_cycle(solver, comm, block_fine, FDstencil, result, max_level, level)
      type(SolverParamsType), intent(in) :: solver
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: block_fine
      type(FDstencil_type), intent(inout) :: FDstencil
      type(ResultType), intent(out) :: result
      integer, intent(in) :: max_level, level

      type(block_type) :: block_coarse

      write(0,*) "Multigrid level: ", level

      if(level == max_level) then
         call iterative_solver(comm, block_fine, FDstencil, solver_params, result)

      else ! If max level is not reached we restrict the residual to a coarser grid

         ! Pre-smoothing
         call iterative_solver(comm, block_fine, FDstencil, solver, result)

         ! Create a coarser grid
         call create_block_type(ndims, 1, 1, domain_begin, domain_end, block_fine%grid_size/2, comm, &
            ghost_begin, ghost_end, stencil_begin, stencil_end, 1, block_coarse)

         ! Residual errors
         call Poisson_residual_2D(0, 1, block_fine, FDstencil)

         ! Restrict the residual to the coarser grid
         call full_weighing_restriction_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%f_matrix_ptr_2D)

         ! We cycle since we are not at the max level
         call w_cycle(solver, comm, block_coarse, FDstencil, result, max_level, level+1)

         ! Prolongate the solution to the finer grid
         call bilinear_prolongation_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%matrix_ptr_2D)

         ! Error correction
         call apply_correction(block_fine%matrix_ptr, block_fine%residual_matrix_ptr)

         ! Second-smoothing
         call iterative_solver(comm, block_fine, FDstencil, solver, result)

         ! Residual again from the second-smoothed solution
         call Poisson_residual_2D(0, 1, block_fine, FDstencil)

         ! Restrict the residual to the coarser grid again
         call full_weighing_restriction_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%f_matrix_ptr_2D)

         ! Second cycle
         call w_cycle(solver, comm, block_coarse, FDstencil, result, max_level, level+1)

         ! Prolongate the solution to the finer grid
         call bilinear_prolongation_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%matrix_ptr_2D)

         ! Final error correction
         call apply_correction(block_fine%matrix_ptr, block_fine%residual_matrix_ptr)

         ! Post-smoothing
         call iterative_solver(comm, block_fine, FDstencil, solver, result)

         ! Deallocate the coarser grid
         call deallocate_block_type(block_coarse)

      end if

   end subroutine w_cycle

   !> Write the solution at a given time
   subroutine write_solution_2D(data_block, t)
      type(block_type), intent(inout) :: data_block
      real, intent(in) :: t

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: u_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(data_block, t) &
      !$omp private(ii, jj, local_indices, global_indices, point, u_val)
      do jj = data_block%extended_block_begin_c(2)+1, data_block%extended_block_end_c(2)
         do ii = data_block%extended_block_begin_c(1)+1, data_block%extended_block_end_c(1)
            local_indices = [ii,jj]
            global_indices = data_block%extended_global_begin_c + local_indices

            call calculate_point(data_block%ndims, -data_block%global_grid_begin, global_indices, &
               data_block%domain_begin, data_block%grid_dx, point)

            call u_analytical_poisson(data_block%ndims, point, t, u_val)

            data_block%matrix_ptr_2D(local_indices(1),local_indices(2)) = u_val

         end do
      end do

      !$omp end parallel do

   end subroutine write_solution_2D

   !> Calculate the rhs of the Poisson equation for a given time
   subroutine write_rhs_2D(data_block, t)
      type(block_type), intent(inout) :: data_block
      real, intent(in) :: t

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: f_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(data_block, t) &
      !$omp private(ii, jj, local_indices, global_indices, point, f_val)
      do jj = data_block%extended_block_begin_c(2)+1, data_block%extended_block_end_c(2)
         do ii = data_block%extended_block_begin_c(1)+1, data_block%extended_block_end_c(1)
            local_indices = [ii,jj]
            global_indices = data_block%extended_global_begin_c + local_indices

            call calculate_point(data_block%ndims, -data_block%global_grid_begin, global_indices, &
               data_block%domain_begin, data_block%grid_dx, point)

            call f_analytical_poisson(ndims, point, t, f_val)

            data_block%f_matrix_ptr_2D(local_indices(1),local_indices(2)) = f_val

         end do
      end do

      !$omp end parallel do

   end subroutine write_rhs_2D

   !> Write the Dirichlet boundary conditions to the solution
   subroutine write_dirichlet_bc_2D(comm_params, data_block)
      type(comm_type), intent(in) :: comm_params
      type(block_type), intent(inout) :: data_block

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: u_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(data_block, current_time) &
      !$omp private(ii, jj, local_indices, global_indices, point, u_val)
      do jj = data_block%extended_block_begin_c(2)+1, data_block%extended_block_end_c(2)
         do ii = data_block%extended_block_begin_c(1)+1, data_block%extended_block_end_c(1)
            local_indices = [ii,jj]
            global_indices = data_block%extended_global_begin_c + local_indices

            call calculate_point(data_block%ndims, -data_block%global_grid_begin, global_indices, &
               data_block%domain_begin, data_block%grid_dx, point)

            call u_analytical_poisson(data_block%ndims, point, current_time, u_val)

            ! At the boundaries of the domain
            if(ii == 1 .or. ii == data_block%extended_block_dims(2) .or. &
               jj == 1 .or. jj == data_block%extended_block_dims(1)) then
               data_block%matrix_ptr_2D(local_indices(1),local_indices(2)) = u_val
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_dirichlet_bc_2D

   subroutine iterative_solver(comm_params, data_block, FDstencil_params, solver_params, result)
      type(comm_type), intent(in) :: comm_params
      type(block_type), intent(inout) :: data_block
      type(FDstencil_type), target, intent(inout) :: FDstencil_params
      type(SolverParamsType), intent(in) :: solver_params
      type(ResultType), intent(out) :: result

      integer :: it, converged

      real, dimension(10) :: norm_array

      call calculate_scaled_coefficients(data_block%ndims, data_block%extended_grid_dx, FDstencil_params)

      call write_dirichlet_bc_2D(comm_params, data_block)
      call sendrecv_data_neighbors(comm_params%comm, data_block, data_block%matrix_ptr)

      norm_array = 1e3

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < solver_params%max_iter)
         if(solver_params%solver_type == 1) then
            call Poisson_Jacobi_2D(omega, 0, 1, data_block, FDstencil_params,&
               norm_array(1), norm_array(2), norm_array(3), norm_array(4))
            call swap_pointers(data_block%matrix_ptr, data_block%temp_matrix_ptr)
            call swap_pointers_2D(data_block%matrix_ptr_2D, data_block%temp_matrix_ptr_2D)
         else if(solver_params%solver_type == 2) then
            call Poisson_Gauss_Seidel_RB_2D(omega, 0, 1, data_block, FDstencil_params,&
               norm_array(1), norm_array(2), norm_array(3), norm_array(4))
            !call Poisson_Gauss_Seidel_2D(omega, 0, 1, data_block, FDstencil_params,&
            !   norm_array(1), norm_array(2), norm_array(3), norm_array(4))
         end if

         call write_dirichlet_bc_2D(comm_params, data_block)
         call sendrecv_data_neighbors(comm_params%comm, data_block, data_block%matrix_ptr)

         call check_convergence(comm_params%comm, solver_params%tol, solver_params%divergence_tol, it, norm_array, converged)
         result%converged = converged
         if(converged == -1) then
            !write(*,*) "Convergence failed"
            !exit
         end if

         it = it + 1
         !converged = 0

      end do

      call set_ResultType(converged, it, norm_array(5), norm_array(6), norm_array(9), norm_array(10), result)

   end subroutine iterative_solver

   !> Write the Dirichlet boundary conditions to the matrix
   subroutine write_matrix_dirichlet_bc_2D(data_block)
      type(block_type), intent(inout) :: data_block

      integer :: ii, jj, diag
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: u_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(data_block, current_time) &
      !$omp private(ii, jj, diag, local_indices, global_indices, point, u_val)
      do jj = data_block%extended_block_begin_c(2)+1, data_block%extended_block_end_c(2)
         do ii = data_block%extended_block_begin_c(1)+1, data_block%extended_block_end_c(1)
            local_indices = [ii,jj]
            global_indices = data_block%extended_global_begin_c + local_indices

            call calculate_point(data_block%ndims, -data_block%global_grid_begin, global_indices, &
               data_block%domain_begin, data_block%grid_dx, point)

            call u_analytical_poisson(data_block%ndims, point, current_time, u_val)

            call IDX_XD(data_block%ndims, data_block%extended_block_dims, local_indices, diag)

            ! At the boundaries of the domain
            if(ii == 1 .or. ii == data_block%extended_block_dims(2) .or. &
               jj == 1 .or. jj == data_block%extended_block_dims(1)) then
               data_block%direct_solver_matrix_ptr_2D(diag,:) = 0.0
               data_block%direct_solver_matrix_ptr_2D(diag,diag) = 1.0
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_matrix_dirichlet_bc_2D

   !> Write the Dirichlet boundary conditions to the rhs for a given time
   subroutine write_rhs_dirichlet_bc_2D(data_block)
      type(block_type), intent(inout) :: data_block

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: u_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(data_block, current_time) &
      !$omp private(ii, jj, local_indices, global_indices, point, u_val)
      do jj = data_block%extended_block_begin_c(2)+1, data_block%extended_block_end_c(2)
         do ii = data_block%extended_block_begin_c(1)+1, data_block%extended_block_end_c(1)
            local_indices = [ii,jj]
            global_indices = data_block%extended_global_begin_c + local_indices

            call calculate_point(data_block%ndims, -data_block%global_grid_begin, global_indices, &
               data_block%domain_begin, data_block%grid_dx, point)

            call u_analytical_poisson(data_block%ndims, point, current_time, u_val)

            ! At the boundaries of the domain
            if(ii == 1 .or. ii == data_block%extended_block_dims(2) .or. &
               jj == 1 .or. jj == data_block%extended_block_dims(1)) then

               data_block%f_matrix_ptr_2D(local_indices(1),local_indices(2)) = u_val
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_rhs_dirichlet_bc_2D

   !!!! POISSON EQUATION FUNCTIONS !!!!

   !> Analytical solution to the Poisson equation
   pure subroutine u_analytical_poisson(ndims, point, t, func_val)
      integer, intent(in) :: ndims
      real, dimension(:), intent(in) :: point
      real, intent(in) :: t
      real, intent(out) :: func_val

      real :: x, y

      x = point(1)
      y = point(2)

      !func_val = product(sin(pi*point))*cos(t)
      !func_val = (sin(pi*x)*sin(2.0*pi*y) + x*sin(pi*y))*cos(t)
      func_val = exp(t)*sin(x)*sin(y)

   end subroutine u_analytical_poisson

   !> Analytical rhs of the Poisson equation u_xx + u_yy = f
   pure subroutine f_analytical_poisson(ndims, point, t, func_val)
      integer, intent(in) :: ndims
      real, dimension(:), intent(in) :: point
      real, intent(in) :: t
      real, intent(out) :: func_val

      real :: x, y

      x = point(1)
      y = point(2)

      !func_val = -ndims*pi*pi*product(sin(pi*point))*cos(t)
      !func_val = -pi*pi*(x + 10.0*sin(pi*x)*cos(pi*y))*sin(pi*y)*cos(t)
      func_val = -2.0*exp(t)*sin(x)*sin(y)

   end subroutine f_analytical_poisson

end module test_Poisson_module



