module poisson_module
   use constants_module, only : pi, MASTER_RANK, filename_txt, filename_dat
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper

   use rank_module, only: rank_type, create_rank_type, deallocate_rank_type, print_rank_type, write_rank_type_blocks_to_file
   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type, print_cartesian_grid
   use FD_module, only: FDstencil_type, create_finite_difference_stencil, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type
   use functions_module, only: FunctionPair, set_function_pointers

   use utility_functions_module, only: open_txt_file, close_txt_file
   use initialization_module, only: write_initial_condition_and_boundary
   use solver_module, only: run_solver
   implicit none

   private
   public :: Poission_1D_analytical, Poission_2D_analytical, Poission_3D_analytical

contains

   !> The analytical solution of the 2D Poisson equation
   subroutine Poission_analytical_wrapper(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      num_derivatives, derivatives, derivatives_sign, alphas, betas)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(:), intent(in) :: grid_size, processor_dims
      real, dimension(:), intent(in) :: domain_begin, domain_end

      integer, intent(in) :: num_derivatives
      integer, dimension(:), intent(in) :: derivatives, alphas, betas
      real, dimension(:), intent(in) :: derivatives_sign

      real, dimension(ndims) :: dx
      integer :: iounit
      real, dimension(8) :: result_array_with_timings

      type(rank_type) :: rank_params
      type(comm_type) :: comm_params
      type(FDstencil_type) :: FDstencil_params
      type(block_type) :: block_params
      type(FunctionPair) :: funcs_params

      ! Setup rank parameters
      call create_rank_type(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, rank_params)

      call set_function_pointers(initial_poisson, boundary_poisson, f_analytical_poisson, u_analytical_poisson, funcs_params)

      call create_cart_comm_type(rank_params%ndims, rank_params%grid_size, rank_params%processor_dim, &
         rank_params%rank, comm_params)

      dx = abs((rank_params%domain_end - rank_params%domain_begin)) / (rank_params%grid_size - 1)

      call create_finite_difference_stencil(rank_params%ndims, num_derivatives, derivatives, &
         derivatives_sign, dx, alphas, betas, FDstencil_params)

      call create_block_type(rank_params%ndims, comm_params, block_params)

      ! Initialize the block
      call write_initial_condition_and_boundary(rank_params%ndims, rank_params%domain_begin, rank_params%domain_end, &
         rank_params%grid_size, block_params%begin, block_params%size, block_params%matrix, &
         FDstencil_params%dx, funcs_params%initial_condition_func, funcs_params%boundary_condition_func)

      ! Time the program
      result_array_with_timings(5) = MPI_WTIME()

      ! Run the solver
      result_array_with_timings(1:4) = run_solver(rank_params, comm_params, block_params, FDstencil_params, funcs_params)

      result_array_with_timings(6) = MPI_WTIME()

      result_array_with_timings(7) = result_array_with_timings(6) - result_array_with_timings(5)

      call all_reduce_mpi_wrapper(result_array_with_timings(7), result_array_with_timings(8), 1, &
         int(MPI_DOUBLE_PRECISION,kind=8), int(MPI_SUM,kind=8), comm_params%comm)

      ! Write out the result and timings from the master rank
      if(rank_params%rank == MASTER_RANK) then

         write(*,"(A, E10.3)") "Glob_norm: ", result_array_with_timings(1)
         write(*,"(A, E10.3)") "Rel_norm: ", result_array_with_timings(2)
         write(*,"(A, F10.1)") "Converged: ", result_array_with_timings(3)
         write(*,"(A, F10.1)") "Iterations: ", result_array_with_timings(4)
         write(*,"(A, F10.3, A)") "Total wall time / processors: ", result_array_with_timings(8)/world_size, " seconds"

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
   subroutine Poission_1D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end)
      integer, intent(in) :: rank, world_size
      integer, dimension(1), intent(in) :: grid_size, processor_dims
      real, dimension(1), intent(in) :: domain_begin, domain_end

      integer, parameter :: num_derivatives = 1
      integer, dimension(num_derivatives), parameter :: derivatives = [2]
      real, dimension(num_derivatives), parameter :: derivatives_sign = [1]
      integer, dimension(1), parameter :: alphas = [1], betas = [1]

      call Poission_analytical_wrapper(1, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, derivatives_sign, alphas, betas)

   end subroutine Poission_1D_analytical

   !> The analytical solution of the 2D Poisson equation
   subroutine Poission_2D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end)
      integer, intent(in) :: rank, world_size
      integer, dimension(2), intent(in) :: grid_size, processor_dims
      real, dimension(2), intent(in) :: domain_begin, domain_end

      integer, parameter :: num_derivatives = 2
      integer, dimension(2*num_derivatives), parameter :: derivatives = [2,0,0,2]
      real, dimension(num_derivatives), parameter :: derivatives_sign = [1,1]
      integer, dimension(2), parameter :: alphas = [1,1], betas = [1,1]

      call Poission_analytical_wrapper(2, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, derivatives_sign, alphas, betas)

   end subroutine Poission_2D_analytical

   !> The analytical solution of the 3D Poisson equation
   subroutine Poission_3D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end)
      integer, intent(in) :: rank, world_size
      integer, dimension(3), intent(in) :: grid_size, processor_dims
      real, dimension(3), intent(in) :: domain_begin, domain_end

      integer, parameter :: num_derivatives = 3
      integer, dimension(3*num_derivatives), parameter :: derivatives = [2,0,0,0,2,0,0,0,2]
      real, dimension(num_derivatives), parameter :: derivatives_sign = [1,1,1]
      integer, dimension(3), parameter :: alphas = [1,1,1], betas = [1,1,1]

      call Poission_analytical_wrapper(3, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, derivatives_sign, alphas, betas)

   end subroutine Poission_3D_analytical

   !!!! POISSON EQUATION FUNCTIONS !!!!

   subroutine initial_poisson(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
      integer, intent(in) :: ndims
      real, dimension(ndims), intent(in) :: domain_start, domain_end, dx, point
      integer, dimension(ndims), intent(in) :: domain_size, domain_indices

      real, intent(inout) :: func_val

      func_val = 20.0

   end subroutine initial_poisson

   subroutine boundary_poisson(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
      integer, intent(in) :: ndims
      real, dimension(ndims), intent(in) :: domain_start, domain_end, dx, point
      integer, dimension(ndims), intent(in) :: domain_size, domain_indices

      real, intent(inout) :: func_val

      integer :: i
      logical :: at_boundary
      at_boundary = .false.

      do i = 1, ndims
         if (domain_indices(i) == 1 .or. domain_indices(i) == domain_size(i)) then
            at_boundary = .true.
            exit
         end if
      end do

      if (at_boundary) then
         call u_analytical_poisson(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
      end if

   end subroutine boundary_poisson

   subroutine u_analytical_poisson(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
      integer, intent(in) :: ndims
      real, dimension(ndims), intent(in) :: domain_start, domain_end, dx, point
      integer, dimension(ndims), intent(in) :: domain_size, domain_indices

      real, intent(inout) :: func_val

      func_val = product(sin(pi*point))

   end subroutine u_analytical_poisson

   subroutine f_analytical_poisson(ndims, domain_start, domain_end, domain_size, domain_indices, dx, point, func_val)
      integer, intent(in) :: ndims
      real, dimension(ndims), intent(in) :: domain_start, domain_end, dx, point
      integer, dimension(ndims), intent(in) :: domain_size, domain_indices

      real, intent(inout) :: func_val

      func_val = -ndims*pi*pi*product(sin(pi*point))

   end subroutine f_analytical_poisson

end module Poisson_module


