module FD_test_module
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

   public :: FD_test_main

contains

   !> x and y dimensions are switched due to FORTRAN column-major ordering. It works if you switch the stencil order at around line 67.
   subroutine FD_test_main(rank, world_size)
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

      grid_size = 8
      processor_dims = 1
      domain_begin = 0
      domain_end = 1
      stencil_sizes = 3

      ! Set the solver parameters tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(1.0, 1.0, 1, 2, solver_params)

      call FD_test(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         stencil_sizes, solver_params)

   end subroutine FD_test_main

   !> The analytical solution of the 2D Navier-Stokes equation
   subroutine FD_test(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      stencil_sizes, solver_params)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(2), intent(in) :: grid_size, processor_dims, stencil_sizes
      real, dimension(2), intent(in) :: domain_begin, domain_end
      type(SolverParamsType), intent(in) :: solver_params

      integer, parameter :: num_derivatives = 4
      integer, dimension(2*num_derivatives), parameter :: derivatives = [0,1,1,0,0,2,2,0]

      call FD_test_test(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, stencil_sizes, solver_params)

   end subroutine FD_test

   !> Solve the 2D Navier-Stokes equation
   subroutine FD_test_test(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
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
      type(FDstencil_type) :: FDstencil_params
      type(block_type) :: test_block_params
      type(FunctionPtrType) :: u_test_2D_func

      u_test_2D_func%output_size = 1
      u_test_2D_func%func => u_test_2D

      call create_cart_comm_type(ndims, processor_dims, rank, world_size, comm_params)

      !call sleeper_function(1)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         grid_size * 0, grid_size * 0, grid_size * 0, grid_size * 0, stencil_sizes/2, stencil_sizes/2, test_block_params)

      call write_function_to_block(test_block_params%ndims, 1, test_block_params%domain_begin, test_block_params%domain_end, &
         test_block_params%extended_grid_size, test_block_params%global_begin_c+1, test_block_params%global_dims, &
         test_block_params%matrix, test_block_params%extended_grid_dx, u_test_2D_func)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      !! Scale the coefficients by dx
      call calculate_scaled_coefficients(ndims, test_block_params%extended_grid_dx, FDstencil_params)

      call test_fd_method(test_block_params, FDstencil_params)

      ! Write out our setup to a file. Just for debugging purposes
      call open_txt_file("output/output_from_", rank, iounit)

      call print_cart_comm_type(comm_params, iounit)

      call print_finite_difference_stencil(ndims, FDstencil_params, iounit)

      call print_block_type(ndims, test_block_params, iounit)

      call close_txt_file(iounit)

      ! ! Write out system solution to a file
      call write_block_data_to_file(test_block_params%data_layout, "output/system_solution.dat", &
         comm_params%comm, test_block_params%matrix)

      ! Deallocate data

      call deallocate_cart_comm_type(comm_params)

      call deallocate_block_type(test_block_params)

      call deallocate_finite_difference_stencil(FDstencil_params)

   end subroutine FD_test_test


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

      dx = test_block%extended_grid_dx

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

            write(*,*) " u_x Diff: ", abs(u_x - u_x_test(1)), " u_y Diff: ", abs(u_y - u_y_test(1)), &
               " u_xx Diff: ", abs(u_xx - u_xx_test(1)), " u_yy Diff: ", abs(u_yy - u_yy_test(1))

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

end module FD_test_module

