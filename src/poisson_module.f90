module Poisson_module
   use constants_module, only : pi, MASTER_RANK, filename_txt, filename_dat
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper

   use rank_module, only: rank_type, create_rank_type, deallocate_rank_type, print_rank_type, write_rank_type_blocks_to_file
   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type, print_cartesian_grid
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, apply_FDstencil
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type
   use functions_module, only: FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD
   use initialization_module, only: write_initial_condition_and_boundary

   use solver_module, only: SolverPtrType, set_SystemSolver_pointer, SolverParamsType, set_SolverParamsType, &
      ResultType, run_solver
   implicit none

   private

   public :: Poission_1D_analytical, Poission_2D_analytical, Poission_3D_analytical

contains

   !> The analytical solution of the 2D Poisson equation
   subroutine Poission_analytical_wrapper(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      num_derivatives, derivatives, stencil_sizes)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(:), intent(in) :: grid_size, processor_dims
      real, dimension(:), intent(in) :: domain_begin, domain_end

      integer, intent(in) :: num_derivatives
      integer, dimension(:), intent(in) :: derivatives, stencil_sizes

      integer, dimension(ndims) ::  begin, end

      real, dimension(ndims) :: dx
      integer :: iounit
      real, dimension(4) :: result_array_with_timings

      type(rank_type) :: rank_params
      type(comm_type) :: comm_params
      type(FDstencil_type) :: FDstencil_params
      type(block_type) :: block_params
      type(FunctionPair) :: funcs_params
      type(SolverPtrType) :: SystemSolver
      type(SolverParamsType) :: solver_params
      type(ResultType) :: result

      call set_SystemSolver_pointer(Poisson_solve_system, SystemSolver)

      ! Set the solver parameters tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(1e-6, 10000, 2, solver_params)

      ! Setup rank parameters
      call create_rank_type(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, rank_params)

      call set_function_pointers(initial_poisson, boundary_poisson, f_analytical_poisson, df_analytical_poisson, &
         u_analytical_poisson, funcs_params)

      call create_cart_comm_type(rank_params%ndims, rank_params%grid_size, rank_params%processor_dim, &
         rank_params%rank, comm_params)

      dx = abs((rank_params%domain_end - rank_params%domain_begin)) / (rank_params%grid_size - 1)

      call create_finite_difference_stencils(rank_params%ndims, num_derivatives, derivatives, &
         dx, stencil_sizes, FDstencil_params)

      call create_block_type(rank_params%ndims, comm_params, block_params)

      ! Initialize the block
      call write_initial_condition_and_boundary(rank_params%ndims, rank_params%domain_begin, rank_params%domain_end, &
         rank_params%grid_size, block_params%begin, block_params%size, block_params%matrix, &
         FDstencil_params%dx, funcs_params)

      ! Set the begin and end for the solver
      begin = 2
      end = block_params%size - 1

      ! Time the program
      result_array_with_timings(1) = MPI_WTIME()

      ! Run the solver
      call run_solver(rank_params, comm_params, block_params, FDstencil_params, funcs_params, SystemSolver, solver_params, &
         begin, end, result)

      result_array_with_timings(2) = MPI_WTIME()

      result_array_with_timings(3) = result_array_with_timings(2) - result_array_with_timings(1)

      call all_reduce_mpi_wrapper(result_array_with_timings(3), result_array_with_timings(4), 1, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_params%comm)

      ! Write out the result and timings from the master rank
      if(rank_params%rank == MASTER_RANK) then
         write(*,"(A, I10.1)") "Converged: ", result%converged
         write(*,"(A, I10.1)") "Iterations: ", result%iterations
         write(*,"(A, E10.3)") "Glob_norm: ", result%global_norm
         write(*,"(A, E10.3)") "Rel_norm: ", result%relative_norm
         write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(4)/world_size, " seconds"

         ! Write out the cartesian grid from the master rank
         call print_cartesian_grid(rank_params%ndims, comm_params%comm, rank_params%world_size, rank_params%processor_dim,&
            filename_txt)
      end if

      ! Write out our setup to a file. Just for debugging purposes
      call open_txt_file(filename_txt, rank_params%rank, iounit)

      call print_rank_type(rank_params, iounit)

      call print_cart_comm_type(rank_params%ndims, comm_params, iounit)

      call print_finite_difference_stencil(rank_params%ndims, FDstencil_params, iounit)

      call print_block_type(rank_params%ndims, block_params, iounit)

      call close_txt_file(iounit)

      ! Write out system solution to a file
      call write_rank_type_blocks_to_file(rank_params, comm_params, block_params, filename_dat)

      ! Deallocate data
      call deallocate_rank_type(rank_params)

      call deallocate_cart_comm_type(comm_params)

      call deallocate_finite_difference_stencil(FDstencil_params)

      call deallocate_block_type(block_params)

   end subroutine Poission_analytical_wrapper

   !> The analytical solution of the 1D Poisson equation
   subroutine Poission_1D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, stencil_sizes)
      integer, intent(in) :: rank, world_size
      integer, dimension(1), intent(in) :: grid_size, processor_dims, stencil_sizes
      real, dimension(1), intent(in) :: domain_begin, domain_end

      integer, parameter :: num_derivatives = 1
      integer, dimension(num_derivatives), parameter :: derivatives = [2]

      call Poission_analytical_wrapper(1, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, stencil_sizes)

   end subroutine Poission_1D_analytical

   !> The analytical solution of the 2D Poisson equation
   subroutine Poission_2D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, stencil_sizes)
      integer, intent(in) :: rank, world_size
      integer, dimension(2), intent(in) :: grid_size, processor_dims, stencil_sizes
      real, dimension(2), intent(in) :: domain_begin, domain_end

      integer, parameter :: num_derivatives = 2
      integer, dimension(2*num_derivatives), parameter :: derivatives = [2,0,0,2]

      call Poission_analytical_wrapper(2, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, stencil_sizes)

   end subroutine Poission_2D_analytical

   !> The analytical solution of the 3D Poisson equation
   subroutine Poission_3D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, stencil_sizes)
      integer, intent(in) :: rank, world_size
      integer, dimension(3), intent(in) :: grid_size, processor_dims, stencil_sizes
      real, dimension(3), intent(in) :: domain_begin, domain_end

      integer, parameter :: num_derivatives = 3
      integer, dimension(3*num_derivatives), parameter :: derivatives = [2,0,0,0,2,0,0,0,2]

      call Poission_analytical_wrapper(3, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, stencil_sizes)

   end subroutine Poission_3D_analytical

   !> Solve the Poisson system
   pure subroutine Poisson_solve_system(ndims, stencil_size, alphas, num_derivatives, stencils, combined_stencils, &
      dims, index, matrix, f_value, F, J)
      integer, intent(in) :: ndims, num_derivatives
      integer, dimension(ndims), intent(in) :: stencil_size, alphas, dims, index
      real, dimension(product(stencil_size)*num_derivatives), intent(in) :: stencils
      real, dimension(product(dims)), intent(in) :: matrix
      real, intent(in) :: f_value

      real, dimension(product(stencil_size)), intent(inout) :: combined_stencils
      real, intent(inout) :: F, J

      integer :: global_index, stencil_begin, stencil_end, num_stencil_elements, center_coefficient_global_index

      combined_stencils = 0.0

      ! Get the number of stencil elements
      num_stencil_elements = product(stencil_size)

      ! Combined stencil into one array since it is linear.
      do global_index = 1, num_derivatives
         stencil_begin = (global_index - 1) * num_stencil_elements + 1
         stencil_end = global_index * num_stencil_elements

         combined_stencils = combined_stencils + stencils(stencil_begin:stencil_end)
      end do

      ! Evaluate the stencil
      call apply_FDstencil(ndims, stencil_size, alphas, combined_stencils, dims, index, matrix, F)
      F = F - f_value

      ! Evaluate the Jacobian. This is simply finding the center coefficient
      call IDX_XD(ndims, stencil_size, alphas + 1, center_coefficient_global_index)
      J = combined_stencils(center_coefficient_global_index)

   end subroutine Poisson_solve_system

   !!!! POISSON EQUATION FUNCTIONS !!!!

   pure subroutine initial_poisson(ndims, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, intent(inout) :: func_val

      func_val = 20.0

   end subroutine initial_poisson

   pure subroutine boundary_poisson(ndims, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, intent(inout) :: func_val

      integer, dimension(ndims) :: block_index
      logical :: at_boundary

      block_index = global_begin_indices + local_indices - 1

      at_boundary = any(block_index == 1 .or. block_index == global_domain_size)

      if (at_boundary) then
         call u_analytical_poisson(ndims, global_begin_indices, local_indices, &
            global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      end if

   end subroutine boundary_poisson

   pure subroutine f_analytical_poisson(ndims, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, intent(inout) :: func_val

      real, dimension(ndims) :: point

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      func_val = -ndims*pi*pi*product(sin(pi*point))

      !func_val = 6.0*point(1) - 4 - 9*pi*pi*sin(3*pi*point(1))

   end subroutine f_analytical_poisson

   pure subroutine df_analytical_poisson(ndims, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, intent(inout) :: func_val

      real, dimension(ndims) :: point

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      func_val = -ndims*pi*pi*pi*product(cos(pi*point))

   end subroutine df_analytical_poisson

   pure subroutine u_analytical_poisson(ndims, global_begin_indices, local_indices, &
      global_domain_begin, global_domain_end, global_domain_size, dx, func_val)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: global_begin_indices, local_indices, global_domain_size
      real, dimension(ndims), intent(in) :: global_domain_begin, global_domain_end, dx

      real, intent(inout) :: func_val

      real, dimension(ndims) :: point

      call calculate_point(ndims, global_begin_indices, local_indices, global_domain_begin, dx, point)

      func_val = product(sin(pi*point))

      !func_val = point(1)*point(1)*point(1) - 2.0*point(1)*point(1) + sin(3*pi*point(1))

   end subroutine u_analytical_poisson

end module Poisson_module



