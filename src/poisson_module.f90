module Poisson_module
   use constants_module, only : pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, apply_FDstencil, update_value_from_stencil, &
      calculate_scaled_coefficients, apply_FDstencil_2D, update_value_from_stencil_2D, set_2D_matrix_coefficients, &
      get_coefficients_wrapper
   use block_module, only: block_type, create_block_type, deallocate_block_type, sendrecv_data_neighbors, print_block_type
   use functions_module, only: FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers, &
      reshape_real_1D_to_2D, print_real_2D_array
   use initialization_module, only: write_initial_condition_and_boundary

   use solver_module, only: SolverPtrType, set_SystemSolver_pointer, SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, choose_iterative_solver, check_convergence
   use multigrid_module, only: full_weighing_restriction_2D, bilinear_prolongation_2D

   implicit none

   private

   !> Number of dimensions and number of derivatives
   integer, parameter :: ndims = 2, num_derivatives = 2
   integer, dimension(ndims*num_derivatives), parameter :: derivatives = [2,0,0,2] ! dxx, dyy

   !> Grid parameters
   integer, dimension(ndims), parameter :: grid_size = [64,64], processor_dims = [1,1]
   logical, dimension(ndims), parameter :: periods = [.false., .false.]
   logical, parameter :: reorder = .true.
   real, dimension(ndims), parameter :: domain_begin = [0,0], domain_end = [1,1]
   integer, dimension(ndims), parameter :: stencil_sizes = 3
   integer, dimension(ndims), parameter :: ghost_begin = [0,0], ghost_end = [0,0]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   !> Solver parameters
   real, parameter :: tol = 1e-36, div_tol = 1e-1, omega = 1.0
   integer, parameter :: max_iter = 30000*(1.0 + 1.0 - omega), multigrid_max_level = 2
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

      integer :: iounit
      real, dimension(4) :: result_array_with_timings

      ! Set the solver parameters tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(tol, div_tol, max_iter, 2, solver_params)

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin, ghost_end, stencil_begin, stencil_end, block_params)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      call calculate_rhs_2D(block_params)
      call write_dirichlet_bc_2D(comm_params, block_params)

      !call sendrecv_data_neighbors(comm_params%comm, block_params, block_params%matrix_ptr)
      !call sendrecv_data_neighbors(comm_params%comm, block_params, block_params%f_matrix_ptr)
      !call sendrecv_data_neighbors(comm_params%comm, block_params, block_params%residual_matrix_ptr)

      !call sleeper_function(1)

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      ! Call the multigrid solver
      call w_cycle(solver_params, comm_params, block_params, FDstencil_params, result, multigrid_max_level, 1)

      !call direct_solver_test(block_params, FDstencil_params)

      call calculate_residual(comm_params, block_params, FDstencil_params)

      result_array_with_timings(2) = MPI_WTIME()

      result_array_with_timings(3) = result_array_with_timings(2) - result_array_with_timings(1)

      call all_reduce_mpi_wrapper(result_array_with_timings(3), result_array_with_timings(4), 1, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_params%comm)

      ! Write out the result and timings from the master rank
      if(rank == MASTER_RANK) then
         call print_resultType(result)
         write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(4)/world_size, " seconds"
      end if

      !call sendrecv_data_neighbors(comm_params%comm, block_params, block_params%matrix_ptr)
      !call sendrecv_data_neighbors(comm_params%comm, block_params, block_params%f_matrix_ptr)
      !call sendrecv_data_neighbors(comm_params%comm, block_params, block_params%residual_matrix_ptr)

      ! Write out our setup to a file. Just for debugging purposes
      call open_txt_file("output/output_from_", rank, iounit)

      !call print_cart_comm_type(comm_params, iounit)

      !call print_finite_difference_stencil(FDstencil_params, iounit)

      call print_block_type(block_params, iounit)

      call close_txt_file(iounit)

      ! Write out system solution to a file
      call write_block_data_to_file(block_params%data_layout, "output/system_solution.dat", &
         comm_params%comm, block_params%matrix)
      call write_block_data_to_file(block_params%data_layout, "output/system_rhs.dat", &
         comm_params%comm, block_params%f_matrix)
      call write_block_data_to_file(block_params%data_layout, "output/system_residual.dat", &
         comm_params%comm, block_params%residual_matrix)

      ! Deallocate data

      call deallocate_cart_comm_type(comm_params)

      call deallocate_block_type(block_params)

      call deallocate_finite_difference_stencil(FDstencil_params)

   end subroutine Poisson_main

   recursive subroutine w_cycle(solver, comm, block_fine, FDstencil, result, max_level, level)
      type(SolverParamsType), intent(in) :: solver
      type(comm_type), intent(in) :: comm
      type(block_type), intent(inout) :: block_fine
      type(FDstencil_type), intent(inout) :: FDstencil
      type(ResultType), intent(out) :: result
      integer, intent(in) :: max_level, level

      type(block_type) :: block_coarse

      write(0,*) "Multigrid level: ", level

      ! Pre-smoothing
      call GS_Method(comm, block_fine, FDstencil, solver, result)

      ! If max level is not reached we restrict the residual to a coarser grid
      if(level /= max_level) then

         ! Create a coarser grid
         call create_block_type(ndims, 1, 1, domain_begin, domain_end, block_fine%grid_size/2, comm, &
            ghost_begin, ghost_end, stencil_begin, stencil_end, block_coarse)

         ! Compute residual errors
         call calculate_residual(comm, block_fine, FDstencil)

         ! Restrict the residual to the coarser grid
         call full_weighing_restriction_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%f_matrix_ptr_2D)

         ! We cycle since we are not at the max level
         call w_cycle(solver, comm, block_coarse, FDstencil, result, max_level, level+1)

         ! Prolongate the solution to the finer grid
         call bilinear_prolongation_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%matrix_ptr_2D)

         ! Error correction
         block_fine%matrix_ptr_2D = block_fine%matrix_ptr_2D + block_fine%residual_matrix_ptr_2D

         ! Second-smoothing
         call GS_Method(comm, block_fine, FDstencil, solver, result)

         ! Calculate residual again
         call calculate_residual(comm, block_fine, FDstencil)

         ! Restrict the residual to the coarser grid again
         call full_weighing_restriction_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%f_matrix_ptr_2D)

         ! Second cycle
         call w_cycle(solver, comm, block_coarse, FDstencil, result, max_level, level+1)

         ! Prolongate the solution to the finer grid
         call bilinear_prolongation_2D(block_fine%extended_block_dims, block_coarse%extended_block_dims, &
            block_fine%residual_matrix_ptr_2D, block_coarse%matrix_ptr_2D)

         ! Final error correction
         block_fine%matrix_ptr_2D = block_fine%matrix_ptr_2D + block_fine%residual_matrix_ptr_2D

         ! Post-smoothing
         call GS_Method(comm, block_fine, FDstencil, solver, result)

         ! Deallocate the coarser grid
         call deallocate_block_type(block_coarse)

      end if

   end subroutine w_cycle

   ! Calculate the residual of the pressure poisson equation
   subroutine calculate_residual(comm, block_params, FDstencil)
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

      !call sendrecv_data_neighbors(comm%comm, block_params, block_params%matrix_ptr)
      !call sendrecv_data_neighbors(comm%comm, block_params, block_params%f_matrix_ptr)

      do ii = block_params%block_begin_c(2)+2, block_params%block_end_c(2)-1
         do jj = block_params%block_begin_c(1)+2, block_params%block_end_c(1)-1
            p_local_indices = [jj,ii]

            call get_coefficients_wrapper(FDstencil, block_params%block_begin_c+1, block_params%block_end_c, p_local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(1:FDstencil%num_stencil_elements)
            dyy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dyy

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

            call apply_FDstencil_2D(combined_stencils_2D, block_params%matrix_ptr_2D, p_local_indices, alpha, beta, laplacian_p)

            block_params%residual_matrix_ptr_2D(jj,ii) = block_params%f_matrix_ptr_2D(jj,ii) - laplacian_p

         end do
      end do

      !call sendrecv_data_neighbors(comm%comm, block_params, block_params%residual_matrix_ptr)

   end subroutine calculate_residual

   subroutine write_dirichlet_bc_2D(comm_params, block_params)
      type(comm_type), intent(in) :: comm_params
      type(block_type), intent(inout) :: block_params

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: u_val

      do ii = block_params%extended_block_begin_c(2)+1, block_params%extended_block_end_c(2)
         do jj = block_params%extended_block_begin_c(1)+1, block_params%extended_block_end_c(1)
            local_indices = [jj,ii]
            global_indices = block_params%extended_global_begin_c + local_indices

            call calculate_point(block_params%ndims, -block_params%global_grid_begin, global_indices, &
               block_params%domain_begin, block_params%grid_dx, point)

            call u_analytical_poisson(block_params%ndims, point, u_val)

            ! At the boundaries of the domain
            if(ii == 1 .or. ii == block_params%extended_block_dims(2) .or. &
               jj == 1 .or. jj == block_params%extended_block_dims(1)) then
               block_params%matrix_ptr_2D(jj,ii) = u_val
            end if

         end do
      end do

      !call sendrecv_data_neighbors(comm_params%comm, block_params, block_params%matrix_ptr)

   end subroutine write_dirichlet_bc_2D

   subroutine GS_Method(comm_params, block_params, FDstencil_params, solver_params, result)
      type(comm_type), intent(in) :: comm_params
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), target, intent(inout) :: FDstencil_params
      type(SolverParamsType), intent(in) :: solver_params
      type(ResultType), intent(out) :: result

      integer :: it, converged

      real, dimension(4) :: norm_array

      call calculate_scaled_coefficients(block_params%ndims, block_params%extended_grid_dx, FDstencil_params)

      call write_dirichlet_bc_2D(comm_params, block_params)

      norm_array = [1e3,1e6,1e9,1e12]

      converged = 0
      it = 0
      do while(converged /= 1 .and. it < solver_params%max_iter)
         call GS_iteration_2D(block_params, FDstencil_params, norm_array(1))

         call write_dirichlet_bc_2D(comm_params, block_params)

         call check_convergence(comm_params%comm, solver_params%tol, solver_params%divergence_tol, &
            1.0/product(block_params%extended_grid_size), norm_array, converged)
         result%converged = converged
         if(converged == -1 .and. it > 0) then
            write(*,*) "Convergence failed"
            exit
         end if

         it = it + 1
         converged = 0

      end do

      result%iterations = it
      result%global_norm = norm_array(3)
      result%relative_norm = norm_array(4)

   end subroutine GS_Method

   subroutine GS_iteration_2D(block_params, FDstencil, local_norm)
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), intent(inout) :: FDstencil
      real, intent(out) :: local_norm

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D
      real :: old_val, f_val, new_val

      do ii = block_params%block_begin_c(2)+2, block_params%block_end_c(2)-1
         do jj = block_params%block_begin_c(1)+2, block_params%block_end_c(1)-1
            local_indices = [jj,ii]

            old_val = block_params%matrix_ptr_2D(jj,ii)
            f_val = block_params%f_matrix_ptr_2D(jj,ii)

            call get_coefficients_wrapper(FDstencil, [1,1], block_params%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dxx => coefficients(1:FDstencil%num_stencil_elements)
            dyy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)

            combined_stencils = dxx + dyy

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, combined_stencils, combined_stencils_2D)

            call update_value_from_stencil_2D(combined_stencils_2D, block_params%matrix_ptr_2D, local_indices, alpha, beta, &
               f_val, new_val)

            new_val = (1.0 - omega) * old_val + omega * new_val

            block_params%matrix_ptr_2D(jj,ii) = new_val

            local_norm = local_norm + (abs(new_val - old_val)**2)

         end do
      end do

   end subroutine GS_iteration_2D

   subroutine calculate_rhs_2D(block_params)
      type(block_type), intent(inout) :: block_params

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: f_val

      do ii = block_params%extended_block_begin_c(2)+1, block_params%extended_block_end_c(2)
         do jj = block_params%extended_block_begin_c(1)+1, block_params%extended_block_end_c(1)
            local_indices = [jj,ii]
            global_indices = block_params%extended_global_begin_c + local_indices

            call calculate_point(block_params%ndims, -block_params%global_grid_begin, global_indices, &
               block_params%domain_begin, block_params%grid_dx, point)

            call f_analytical_poisson(ndims, point, f_val)

            block_params%f_matrix_ptr_2D(jj,ii) = f_val

         end do
      end do

   end subroutine calculate_rhs_2D

   subroutine direct_solver_test(block_params, FDstencil_params)
      type(block_type), intent(inout) :: block_params
      type(FDstencil_type), intent(inout) :: FDstencil_params

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, alpha, beta
      real, contiguous, dimension(:), pointer :: coefficients, dxx, dyy
      real, dimension(product(FDstencil_params%stencil_sizes)), target :: combined_stencils
      real, contiguous, dimension(:,:), pointer :: combined_stencils_2D

      integer :: num_equations, num_rhs, lda, ldb
      integer, dimension(size(block_params%direct_solver_matrix_ptr_2D,1)) :: ipiv
      integer :: info

      call calculate_scaled_coefficients(block_params%ndims, block_params%extended_grid_dx, FDstencil_params)

      do ii = block_params%block_begin_c(2)+1, block_params%block_end_c(2)
         do jj = block_params%block_begin_c(1)+1, block_params%block_end_c(1)
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

      ! Calculate the right hand side
      call calculate_rhs_2D(block_params)
      ! Apply the Dirichlet boundary conditions
      call write_direct_dirichlet_bc_2D(block_params)

      num_equations = size(block_params%direct_solver_matrix_ptr_2D,1)
      num_rhs = 1
      lda = size(block_params%direct_solver_matrix_ptr_2D,1)
      ldb = size(block_params%f_matrix_ptr,1)

      info = 0
      call DGESV(num_equations, num_rhs, block_params%direct_solver_matrix_ptr_2D, lda, ipiv, &
         block_params%f_matrix_ptr, ldb, info)
      if (info /= 0) then
         print *, "DGESV reported an error: ", info
         stop
      end if

   end subroutine direct_solver_test

   subroutine write_direct_dirichlet_bc_2D(block_params)
      type(block_type), intent(inout) :: block_params

      integer :: ii, jj, diag
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real :: u_val

      do ii = block_params%block_begin_c(2)+1, block_params%block_end_c(2)
         do jj = block_params%block_begin_c(1)+1, block_params%block_end_c(1)
            local_indices = [jj,ii]
            global_indices = block_params%extended_global_begin_c + local_indices

            call calculate_point(block_params%ndims, -block_params%global_grid_begin, global_indices, &
               block_params%domain_begin, block_params%grid_dx, point)

            call u_analytical_poisson(block_params%ndims, point, u_val)

            diag = local_indices(1) + (local_indices(2) - 1) * block_params%extended_block_dims(1)

            ! At the boundaries of the domain
            if(ii == 1 .or. ii == block_params%extended_block_dims(2) .or. &
               jj == 1 .or. jj == block_params%extended_block_dims(1)) then
               block_params%direct_solver_matrix_ptr_2D(diag,:) = 0.0 ! Diag in the first column?
               block_params%direct_solver_matrix_ptr_2D(diag,diag) = 1.0
               block_params%f_matrix_ptr_2D(jj,ii) = u_val
            end if

         end do
      end do

   end subroutine write_direct_dirichlet_bc_2D

   !!!! POISSON EQUATION FUNCTIONS !!!!

   pure subroutine u_analytical_poisson(ndims, point, func_val)
      integer, intent(in) :: ndims
      real, dimension(:), intent(in) :: point
      real, intent(out) :: func_val

      real :: x, y

      x = point(1)
      y = point(2)

      !func_val = product(sin(pi*point))
      func_val = sin(pi*x)*sin(2.0*pi*y) + x*sin(pi*y)

   end subroutine u_analytical_poisson

   pure subroutine f_analytical_poisson(ndims, point, func_val)
      integer, intent(in) :: ndims
      real, dimension(:), intent(in) :: point
      real, intent(out) :: func_val

      real :: x, y

      x = point(1)
      y = point(2)

      !func_val = -ndims*pi*pi*product(sin(pi*point))
      func_val = -pi*pi*(x + 10.0*sin(pi*x)*cos(pi*y))*sin(pi*y)

   end subroutine f_analytical_poisson

end module Poisson_module



