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
   use functions_module, only: FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers, &
      swap_pointers_2D, reshape_real_1D_to_2D, print_real_2D_array, calculate_dt_from_CFL

   use solver_module, only: SolverPtrType, set_SystemSolver_pointer, SolverParamsType, set_SolverParamsType, &
      ResultType, set_ResultType, print_ResultType, choose_iterative_solver, check_convergence, &
      LU_decomposition, solve_LU_system
   use multigrid_module, only: full_weighing_restriction_2D, bilinear_prolongation_2D, apply_correction_2D

   implicit none

   private

   !> Number of dimensions and number of derivatives
   integer, parameter :: ndims = 2, num_derivatives = 2
   integer, dimension(ndims*num_derivatives), parameter :: derivatives = [2,0,0,2] ! dyy, dxx

   !> Grid parameters
   integer, dimension(ndims), parameter :: grid_size = [32,32], processor_dims = [1,1]
   logical, dimension(ndims), parameter :: periods = [.false., .false.]
   logical, parameter :: reorder = .true.
   real, dimension(ndims), parameter :: domain_begin = [0,0], domain_end = [pi,pi]
   integer, dimension(ndims), parameter :: stencil_sizes = 5
   integer, dimension(ndims), parameter :: ghost_begin = [0,0], ghost_end = [0,0]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   ! Time parameters
   real, parameter :: t_0 = 0.0, t_1 = 1.0
   real :: current_time = t_0
   real :: t_steps, dt

   !> Solver parameters
   integer, parameter :: direct_or_iterative = 0, Jacobi_or_GS = 1
   real, parameter :: tol = (1e-12)**2, div_tol = 1e-1, omega = 0.8
   integer, parameter :: max_iter = 10000 * (2.0 - omega), multigrid_max_level = 2
   type(SolverParamsType) :: solver_params

   public :: Poisson_main

contains

   !> The main Poisson subroutine to solve the system
   subroutine Poisson_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      type(comm_type) :: comm_params
      type(FDstencil_type) :: FDstencil_params
      type(block_type) :: block_params
      type(ResultType) :: result

      integer :: iounit, ii
      real, dimension(4) :: result_array_with_timings

      !call sleeper_function(1)

      ! Set the solver parameters tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(tol, div_tol, max_iter, Jacobi_or_GS, solver_params)

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin, ghost_end, stencil_begin, stencil_end, direct_or_iterative, block_params)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      !call calculate_dt_from_CFL(1.0, [1.0,1.0], block_params%extended_grid_dx, dt)
      !t_steps = (t_end - t_begin) / dt + 1
      t_steps = 1
      dt = (t_1 - t_0) / t_steps

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      if(direct_or_iterative == 1) then
         do ii = 1, t_steps
            current_time = current_time + dt
            call one_timestep(block_params, FDstencil_params, comm_params, result)
            if(rank == MASTER_RANK) then
               call print_ResultType(result)
               print*, "Time step: ", ii, " Time: ", current_time
            end if
         end do

      else
         !> Assemble the matrix A
         call assemble_matrix(block_params, FDstencil_params)

         ! Apply the Dirichlet boundary conditions to the matrix before decomposing it
         call write_matrix_dirichlet_bc(block_params)

         !> Decompose the matrix A into LU
         call LU_decomposition(block_params)

         !> Write the initial condition to the system
         call write_solution_2D(block_params, t_0)

         do ii = 1, t_steps
            current_time = current_time + dt
            call one_timestep(block_params, FDstencil_params, comm_params, result)
            if(rank == MASTER_RANK) then
               print*, "Time step: ", ii, " Time: ", current_time
            end if
         end do

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

      !call print_finite_difference_stencil(FDstencil_params, iounit)

      call print_block_type(block_params, iounit)

      call close_txt_file(iounit)

      ! Write out system solution to a file
      if(direct_or_iterative == 1) then
         call write_block_data_to_file(block_params%data_layout, "output/system_solution.dat", &
            comm_params%comm, block_params%matrix)
         call write_block_data_to_file(block_params%data_layout, "output/system_rhs.dat", &
            comm_params%comm, block_params%f_matrix)
         call write_block_data_to_file(block_params%data_layout, "output/system_residual.dat", &
            comm_params%comm, block_params%residual_matrix)
      else
         call write_block_data_to_file(block_params%data_layout, "output/system_solution.dat", &
            comm_params%comm, block_params%f_matrix)
         call write_block_data_to_file(block_params%data_layout, "output/system_rhs.dat", &
            comm_params%comm, block_params%matrix_ptr)
         call write_block_data_to_file(block_params%data_layout, "output/system_residual.dat", &
            comm_params%comm, block_params%f_matrix)
      end if

      ! Deallocate data

      call deallocate_cart_comm_type(comm_params)

      call deallocate_block_type(block_params)

      call deallocate_finite_difference_stencil(FDstencil_params)

   end subroutine Poisson_main

   subroutine one_timestep(block_params, FDstencil_params, comm_params, result)
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), intent(inout) :: FDstencil_params
      type(comm_type), intent(in) :: comm_params
      type(ResultType), intent(out) :: result

      !> Calculate the rhs of the system
      call write_rhs_2D(block_params, current_time)

      if(direct_or_iterative == 1) then
         !> Apply the Dirichlet boundary conditions to the system solution
         call write_dirichlet_bc_2D(comm_params, block_params)

         ! Call the multigrid solver
         call w_cycle(solver_params, comm_params, block_params, FDstencil_params, result, multigrid_max_level, 1)

         block_params%new_matrix_ptr_2D = block_params%new_matrix_ptr_2D + (dt/-2.0) * block_params%matrix_ptr_2D ! Parallelize this

         call swap_pointers(block_params%matrix_ptr, block_params%new_matrix_ptr)
      else
         !> Apply the Dirichlet boundary conditions to the rhs
         call write_rhs_dirichlet_bc(block_params)

         !> Solve the linear system using the LU decomposition
         call solve_LU_system(block_params)

         block_params%new_matrix_ptr_2D = block_params%new_matrix_ptr_2D + (dt/-2.0) * block_params%f_matrix_ptr_2D
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
         call calculate_residual_2D(comm, block_fine, FDstencil)

         ! Restrict the residual to the coarser grid
         call full_weighing_restriction_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%f_matrix_ptr_2D)

         ! We cycle since we are not at the max level
         call w_cycle(solver, comm, block_coarse, FDstencil, result, max_level, level+1)

         ! Prolongate the solution to the finer grid
         call bilinear_prolongation_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%matrix_ptr_2D)

         ! Error correction
         call apply_correction_2D(block_fine%extended_block_dims,block_fine%matrix_ptr_2D, block_fine%residual_matrix_ptr_2D)

         ! Second-smoothing
         call iterative_solver(comm, block_fine, FDstencil, solver, result)

         ! Residual again from the second-smoothed solution
         call calculate_residual_2D(comm, block_fine, FDstencil)

         ! Restrict the residual to the coarser grid again
         call full_weighing_restriction_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%f_matrix_ptr_2D)

         ! Second cycle
         call w_cycle(solver, comm, block_coarse, FDstencil, result, max_level, level+1)

         ! Prolongate the solution to the finer grid
         call bilinear_prolongation_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%matrix_ptr_2D)

         ! Final error correction
         call apply_correction_2D(block_fine%extended_block_dims,block_fine%matrix_ptr_2D, block_fine%residual_matrix_ptr_2D)

         ! Post-smoothing
         call iterative_solver(comm, block_fine, FDstencil, solver, result)

         ! Deallocate the coarser grid
         call deallocate_block_type(block_coarse)

      end if

   end subroutine w_cycle

   !> Write the initial condition to the system
   subroutine write_solution_2D(block_params, t)
      type(block_type), intent(inout) :: block_params
      real, intent(in) :: t

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: u_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(block_params, t) &
      !$omp private(ii, jj, local_indices, global_indices, point, u_val)
      do ii = block_params%extended_block_begin_c(2)+1, block_params%extended_block_end_c(2)
         do jj = block_params%extended_block_begin_c(1)+1, block_params%extended_block_end_c(1)
            local_indices = [jj,ii]
            global_indices = block_params%extended_global_begin_c + local_indices

            call calculate_point(block_params%ndims, -block_params%global_grid_begin, global_indices, &
               block_params%domain_begin, block_params%grid_dx, point)

            call u_analytical_poisson(block_params%ndims, point, t, u_val)

            block_params%matrix_ptr_2D(jj,ii) = u_val

         end do
      end do

      !$omp end parallel do

   end subroutine write_solution_2D

   !> Calculate the rhs of the Poisson equation for a given time
   subroutine write_rhs_2D(block_params, t)
      type(block_type), intent(inout) :: block_params
      real, intent(in) :: t

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: f_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(block_params, t) &
      !$omp private(ii, jj, local_indices, global_indices, point, f_val)
      do ii = block_params%extended_block_begin_c(2)+1, block_params%extended_block_end_c(2)
         do jj = block_params%extended_block_begin_c(1)+1, block_params%extended_block_end_c(1)
            local_indices = [jj,ii]
            global_indices = block_params%extended_global_begin_c + local_indices

            call calculate_point(block_params%ndims, -block_params%global_grid_begin, global_indices, &
               block_params%domain_begin, block_params%grid_dx, point)

            call f_analytical_poisson(ndims, point, t, f_val)

            block_params%f_matrix_ptr_2D(jj,ii) = f_val

         end do
      end do

      !$omp end parallel do

   end subroutine write_rhs_2D

   !> Calculate the residual of the pressure poisson equation
   subroutine calculate_residual_2D(comm, block_params, FDstencil)
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), intent(inout) :: FDstencil

      integer :: ii, jj
      integer, dimension(ndims) :: p_local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
      real :: laplacian_p

      call calculate_scaled_coefficients(block_params%ndims, block_params%extended_grid_dx, FDstencil)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(block_params, FDstencil) &
      !$omp private(ii, jj, p_local_indices, alpha, beta, coefficients, dxx, dyy, combined_stencils, &
      !$omp combined_stencils_2D, laplacian_p)
      do ii = block_params%block_begin_c(2)+1, block_params%block_end_c(2)
         do jj = block_params%block_begin_c(1)+1, block_params%block_end_c(1)
            p_local_indices = [jj,ii]

            call get_coefficients_wrapper(FDstencil, [1,1], block_params%extended_block_dims, p_local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(1:FDstencil%num_stencil_elements)
            dyy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dyy

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

            call apply_FDstencil_2D(combined_stencils_2D, block_params%matrix_ptr_2D, p_local_indices, alpha, beta, laplacian_p)

            block_params%residual_matrix_ptr_2D(jj,ii) = block_params%f_matrix_ptr_2D(jj,ii) - laplacian_p

         end do
      end do

      !$omp end parallel do

      ! Depends on the boundary conditions this is for Dirichlet. For Neumann we need to calculate the gradient at the point. Parallelize this
      block_params%residual_matrix_ptr_2D(:,1) = 0.0
      block_params%residual_matrix_ptr_2D(:,size(block_params%residual_matrix_ptr_2D,2)) = 0.0
      block_params%residual_matrix_ptr_2D(1,:) = 0.0
      block_params%residual_matrix_ptr_2D(size(block_params%residual_matrix_ptr_2D,1),:) = 0.0

   end subroutine calculate_residual_2D

   !> Write the Dirichlet boundary conditions to the solution
   subroutine write_dirichlet_bc_2D(comm_params, block_params)
      type(comm_type), intent(in) :: comm_params
      type(block_type), intent(inout) :: block_params

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: u_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(block_params, current_time) &
      !$omp private(ii, jj, local_indices, global_indices, point, u_val)
      do ii = block_params%extended_block_begin_c(2)+1, block_params%extended_block_end_c(2)
         do jj = block_params%extended_block_begin_c(1)+1, block_params%extended_block_end_c(1)
            local_indices = [jj,ii]
            global_indices = block_params%extended_global_begin_c + local_indices

            call calculate_point(block_params%ndims, -block_params%global_grid_begin, global_indices, &
               block_params%domain_begin, block_params%grid_dx, point)

            call u_analytical_poisson(block_params%ndims, point, current_time, u_val)

            ! At the boundaries of the domain
            if(ii == 1 .or. ii == block_params%extended_block_dims(2) .or. &
               jj == 1 .or. jj == block_params%extended_block_dims(1)) then
               block_params%matrix_ptr_2D(jj,ii) = u_val
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_dirichlet_bc_2D

   subroutine iterative_solver(comm_params, block_params, FDstencil_params, solver_params, result)
      type(comm_type), intent(in) :: comm_params
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), target, intent(inout) :: FDstencil_params
      type(SolverParamsType), intent(in) :: solver_params
      type(ResultType), intent(out) :: result

      integer :: it, converged

      real, dimension(10) :: norm_array

      call calculate_scaled_coefficients(block_params%ndims, block_params%extended_grid_dx, FDstencil_params)

      call write_dirichlet_bc_2D(comm_params, block_params)

      norm_array = 1e3

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < solver_params%max_iter)
         if(solver_params%solver_type == 1) then
            call Jacobi_iteration_2D(block_params, FDstencil_params, norm_array(1), norm_array(2), norm_array(3), norm_array(4))
            call swap_pointers_2D(block_params%matrix_ptr_2D, block_params%temp_matrix_ptr_2D)
         else if(solver_params%solver_type == 2) then
            call GS_iteration_2D(block_params, FDstencil_params, norm_array(1), norm_array(2), norm_array(3), norm_array(4))
         end if

         call write_dirichlet_bc_2D(comm_params, block_params)

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

   subroutine GS_iteration_2D(block_params, FDstencil, local_u_diff_norm, local_r_diff_norm, local_u_norm, local_r_norm)
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), intent(inout) :: FDstencil
      real, intent(out) :: local_u_diff_norm, local_r_diff_norm, local_u_norm, local_r_norm

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
      real :: f_val, u0_val, u1_val, r0_val, r1_val

      local_u_diff_norm = 0.0
      local_r_diff_norm = 0.0
      local_u_norm = 0.0
      local_r_norm = 0.0

      ! Do red black scheme here
      do ii = block_params%block_begin_c(2)+1, block_params%block_end_c(2)
         do jj = block_params%block_begin_c(1)+1, block_params%block_end_c(1)
            local_indices = [jj,ii]

            f_val = block_params%f_matrix_ptr_2D(jj,ii)
            u0_val = block_params%matrix_ptr_2D(jj,ii)
            r0_val = block_params%residual_matrix_ptr_2D(jj,ii)

            call get_coefficients_wrapper(FDstencil, [1,1], block_params%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(1:FDstencil%num_stencil_elements)
            dyy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dyy

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

            call update_value_from_stencil_2D(combined_stencils_2D, block_params%matrix_ptr_2D, local_indices, alpha, beta, &
               f_val, u1_val, r1_val)

            u1_val = (1.0 - omega) * u0_val + omega * u1_val

            block_params%matrix_ptr_2D(jj,ii) = u1_val
            local_u_diff_norm = local_u_diff_norm + (abs(u1_val - u0_val)**2)
            local_u_norm = local_u_norm + (abs(u1_val)**2)

         end do
      end do

   end subroutine GS_iteration_2D

   subroutine Jacobi_iteration_2D(block_params, FDstencil, local_u_diff_norm, local_r_diff_norm, local_u_norm, local_r_norm)
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), intent(inout) :: FDstencil
      real, intent(out) :: local_u_diff_norm, local_r_diff_norm, local_u_norm, local_r_norm

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
      real :: f_val, u0_val, u1_val, r0_val, r1_val

      local_u_diff_norm = 0.0
      local_r_diff_norm = 0.0
      local_u_norm = 0.0
      local_r_norm = 0.0

      !$omp parallel do collapse(2) default(none) reduction(+:local_u_diff_norm, local_u_norm) &
      !$omp shared(block_params, FDstencil) &
      !$omp private(ii, jj, local_indices, alpha, beta, coefficients, dxx, dyy, combined_stencils, &
      !$omp combined_stencils_2D, f_val, u0_val, u1_val, r0_val, r1_val)
      do ii = block_params%block_begin_c(2)+1, block_params%block_end_c(2)
         do jj = block_params%block_begin_c(1)+1, block_params%block_end_c(1)
            local_indices = [jj,ii]

            f_val = block_params%f_matrix_ptr_2D(jj,ii)
            u0_val = block_params%matrix_ptr_2D(jj,ii)
            r0_val = block_params%residual_matrix_ptr_2D(jj,ii)

            call get_coefficients_wrapper(FDstencil, [1,1], block_params%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(1:FDstencil%num_stencil_elements)
            dyy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dyy

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

            call update_value_from_stencil_2D(combined_stencils_2D, block_params%matrix_ptr_2D, local_indices, alpha, beta, &
               f_val, u1_val, r1_val)

            u1_val = (1.0 - omega) * u0_val + omega * u1_val

            block_params%temp_matrix_ptr_2D(jj,ii) = u1_val
            local_u_diff_norm = local_u_diff_norm + (abs(u1_val - u0_val)**2)
            local_u_norm = local_u_norm + (abs(u1_val)**2)

         end do
      end do

      !$omp end parallel do

   end subroutine Jacobi_iteration_2D

   !> The global assembly of the matrix A. The vector f is already assembled
   subroutine assemble_matrix(block_params, FDstencil_params)
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), intent(inout) :: FDstencil_params

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil_params%stencil_sizes)), target :: combined_stencils

      call calculate_scaled_coefficients(block_params%ndims, block_params%extended_grid_dx, FDstencil_params)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(block_params, FDstencil_params) &
      !$omp private(ii, jj, local_indices, alpha, beta, coefficients, dxx, dyy, combined_stencils)
      do ii = block_params%extended_block_begin_c(2)+1, block_params%extended_block_end_c(2)
         do jj = block_params%extended_block_begin_c(1)+1, block_params%extended_block_end_c(1)
            local_indices = [jj,ii]

            call get_coefficients_wrapper(FDstencil_params, [1,1], block_params%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(1:FDstencil_params%num_stencil_elements)
            dyy => coefficients(FDstencil_params%num_stencil_elements + 1:2 * FDstencil_params%num_stencil_elements)

            combined_stencils = dxx + dyy

            call set_matrix_coefficients(block_params%ndims, FDstencil_params%stencil_sizes, block_params%extended_block_dims, &
               combined_stencils, block_params%direct_solver_matrix_ptr_2D, local_indices, alpha, beta)

         end do
      end do

      !$omp end parallel do

   end subroutine assemble_matrix

   !> Write the Dirichlet boundary conditions to the matrix
   subroutine write_matrix_dirichlet_bc(block_params)
      type(block_type), intent(inout) :: block_params

      integer :: ii, jj, diag
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: u_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(block_params, current_time) &
      !$omp private(ii, jj, diag, local_indices, global_indices, point, u_val)
      do ii = block_params%extended_block_begin_c(2)+1, block_params%extended_block_end_c(2)
         do jj = block_params%extended_block_begin_c(1)+1, block_params%extended_block_end_c(1)
            local_indices = [jj,ii]
            global_indices = block_params%extended_global_begin_c + local_indices

            call calculate_point(block_params%ndims, -block_params%global_grid_begin, global_indices, &
               block_params%domain_begin, block_params%grid_dx, point)

            call u_analytical_poisson(block_params%ndims, point, current_time, u_val)

            call IDX_XD(block_params%ndims, block_params%extended_block_dims, local_indices, diag)

            ! At the boundaries of the domain
            if(ii == 1 .or. ii == block_params%extended_block_dims(2) .or. &
               jj == 1 .or. jj == block_params%extended_block_dims(1)) then
               block_params%direct_solver_matrix_ptr_2D(diag,:) = 0.0
               block_params%direct_solver_matrix_ptr_2D(diag,diag) = 1.0
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_matrix_dirichlet_bc

   !> Write the Dirichlet boundary conditions to the rhs for a given time
   subroutine write_rhs_dirichlet_bc(block_params)
      type(block_type), intent(inout) :: block_params

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: u_val

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(block_params, current_time) &
      !$omp private(ii, jj, local_indices, global_indices, point, u_val)
      do ii = block_params%extended_block_begin_c(2)+1, block_params%extended_block_end_c(2)
         do jj = block_params%extended_block_begin_c(1)+1, block_params%extended_block_end_c(1)
            local_indices = [jj,ii]
            global_indices = block_params%extended_global_begin_c + local_indices

            call calculate_point(block_params%ndims, -block_params%global_grid_begin, global_indices, &
               block_params%domain_begin, block_params%grid_dx, point)

            call u_analytical_poisson(block_params%ndims, point, current_time, u_val)

            ! At the boundaries of the domain
            if(ii == 1 .or. ii == block_params%extended_block_dims(2) .or. &
               jj == 1 .or. jj == block_params%extended_block_dims(1)) then

               block_params%f_matrix_ptr_2D(jj,ii) = u_val
            end if

         end do
      end do

      !$omp end parallel do

   end subroutine write_rhs_dirichlet_bc

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



