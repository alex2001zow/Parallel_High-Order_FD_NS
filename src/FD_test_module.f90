module FD_test_module
   use constants_module, only: pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, apply_FDstencil_2D, calculate_scaled_coefficients, get_coefficients_wrapper
   use block_module, only: block_type, create_block_type, deallocate_block_type, sendrecv_data_neighbors, print_block_type
   use functions_module, only: calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function, swap_pointers, &
      reshape_real_1D_to_2D

   use solver_module, only: SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, check_convergence
   implicit none

   private

   !> Number of dimensions of number of derivatives
   integer, parameter :: ndims = 2, num_derivatives = 5
   integer, dimension(ndims * num_derivatives), parameter :: derivatives = [1,0,0,1,1,1,2,0,0,2] ! dx, dy, dxy, dxx, dyy

   !> Grid parameters
   integer, dimension(ndims), parameter :: grid_size = [128,128]
   integer, dimension(ndims), parameter :: processor_dims = [1,1]
   logical, dimension(ndims), parameter :: periods = [.false., .false.]
   logical, parameter :: reorder = .true.
   real, dimension(ndims), parameter :: domain_begin = [-2,-2], domain_end = [2,2]
   integer, dimension(ndims), parameter :: stencil_sizes = 3
   integer, dimension(ndims), parameter :: ghost_begin = [0,0] , ghost_end = [0,0]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end=stencil_sizes/2

   public :: FD_test_main

contains

   subroutine FD_test_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      integer :: iounit

      type(comm_type) :: comm_params
      type(FDstencil_type) :: FDstencil_params
      type(block_type) :: FD_block

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm_params)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin, ghost_end, stencil_begin, stencil_end, 1, FD_block)

      call write_analytical_function(FD_block)

      call create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_params)

      ! Apply the FD stencil to the block
      call test_fd_method(FD_block, FDstencil_params)

      ! Write out our setup to a file. Just for debugging purposes
      !call open_txt_file("output/output_from_", rank, iounit)

      !call print_cart_comm_type(comm_params, iounit)

      !call print_finite_difference_stencil(FDstencil_params, iounit)

      !call print_block_type(FD_block, iounit)

      !call close_txt_file(iounit)

      ! Write out system solution to a file
      call write_block_data_to_file(FD_block%data_layout, "output/system_solution.dat", &
         comm_params%comm, FD_block%matrix_ptr)

      call write_block_data_to_file(FD_block%data_layout, "output/residual_solution.dat", &
         comm_params%comm, FD_block%f_matrix_ptr)

      ! Deallocate data
      call deallocate_cart_comm_type(comm_params)

      call deallocate_block_type(FD_block)

      call deallocate_finite_difference_stencil(FDstencil_params)

   end subroutine FD_test_main

   subroutine test_fd_method(test_block, FDstencil)
      type(block_type), intent(inout) :: test_block
      type(FDstencil_type), intent(inout) :: FDstencil

      integer :: ii, jj, i, j

      real, contiguous, dimension(:), pointer :: coefficients
      real, contiguous, dimension(:), pointer :: dx, dy, dxy, dxx, dyy
      real, contiguous, dimension(:,:), pointer :: dx_2D, dy_2D, dxy_2D, dxx_2D, dyy_2D
      integer, dimension(ndims) :: stencil_size, alpha, beta, start_dims, local_indices

      real :: u_x, u_y, u_xy, u_xx, u_yy

      call calculate_scaled_coefficients(ndims, test_block%extended_grid_dx, FDstencil)

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(test_block, FDstencil) &
      !$omp private(ii, jj, local_indices, coefficients, dx, dy, dxy, dxx, dyy, dx_2D, dy_2D, dxy_2D, dxx_2D, dyy_2D, &
      !$omp alpha, beta, u_x, u_y, u_xy, u_xx, u_yy)
      do jj = test_block%block_begin_c(2)+1, test_block%block_end_c(2)
         do ii = test_block%block_begin_c(1)+1, test_block%block_end_c(1)
            local_indices = [ii,jj]

            call get_coefficients_wrapper(FDstencil, [1,1], test_block%extended_block_dims, local_indices, &
               alpha, beta, coefficients)

            dx => coefficients(1:FDstencil%num_stencil_elements)
            dy => coefficients(FDstencil%num_stencil_elements + 1:2 * FDstencil%num_stencil_elements)
            dxy => coefficients(2 * FDstencil%num_stencil_elements + 1:3 * FDstencil%num_stencil_elements)
            dxx => coefficients(3 * FDstencil%num_stencil_elements + 1:4 * FDstencil%num_stencil_elements)
            dyy => coefficients(4 * FDstencil%num_stencil_elements + 1:5 * FDstencil%num_stencil_elements)

            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dx, dx_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dy, dy_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dxy, dxy_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dxx, dxx_2D)
            call reshape_real_1D_to_2D(FDstencil%stencil_sizes, dyy, dyy_2D)

            call apply_FDstencil_2D(dx_2D, test_block%matrix_ptr_2D, local_indices, alpha, beta, u_x)
            call apply_FDstencil_2D(dy_2D, test_block%matrix_ptr_2D, local_indices, alpha, beta, u_y)
            call apply_FDstencil_2D(dxy_2D, test_block%matrix_ptr_2D, local_indices, alpha, beta, u_xy)
            call apply_FDstencil_2D(dxx_2D, test_block%matrix_ptr_2D, local_indices, alpha, beta, u_xx)
            call apply_FDstencil_2D(dyy_2D, test_block%matrix_ptr_2D, local_indices, alpha, beta, u_yy)

            test_block%f_matrix_ptr_2D(local_indices(1), local_indices(2)) = u_xy

         end do
      end do

      !$omp end parallel do

   end subroutine test_fd_method

   subroutine write_analytical_function(test_block)
      type (block_type), intent(inout) :: test_block

      integer :: ii, jj
      integer, dimension(ndims) :: local_indices, global_indices
      real, dimension(ndims) :: point
      real ::u, u_x, u_y, u_xy, u_xx, u_yy

      !$omp parallel do collapse(2) default(none) &
      !$omp shared(test_block) &
      !$omp private(ii, jj, local_indices, global_indices, point, u, u_x, u_y, u_xy, u_xx, u_yy)
      do jj = test_block%block_begin_c(2)+1, test_block%block_end_c(2)
         do ii = test_block%block_begin_c(1)+1, test_block%block_end_c(1)
            local_indices = [ii,jj]
            global_indices = test_block%extended_global_begin_c + local_indices

            call calculate_point(test_block%ndims, -test_block%global_grid_begin, global_indices, &
               test_block%domain_begin, test_block%grid_dx, point)

            call analytical_test_function_2D(point, u, u_x, u_y, u_xy, u_xx, u_yy)

            test_block%matrix_ptr_2D(local_indices(1), local_indices(2)) = u

         end do
      end do

      !$omp end parallel do

   end subroutine write_analytical_function

   pure subroutine analytical_test_function_2D(point, u, u_x, u_y, u_xy, u_xx, u_yy)
      real, dimension(2), intent(in) :: point
      real, intent(out) :: u, u_x, u_y, u_xy, u_xx, u_yy

      real :: x, y

      x = point(1)
      y = point(2)

      ! u = sin(x)*cos(y)
      ! u_x = cos(x)*cos(y)
      ! u_y = -sin(x)*sin(y)
      ! u_xy = -sin(x)*cos(y)
      ! u_xx = -sin(x)*cos(y)
      ! u_yy = -sin(x)*cos(y)

      u = 3*x*x + 2*x*y*y + y*y*y + y*y*exp(3*x) - 4*y + sin(2*y)*cos(3*x) + 7
      u_x = 6*x + 3*y*y*exp(3*x) + 2*y*y - 3*sin(2*y)*sin(3*x)
      u_y = 4*x*y + 3*y*y + 2*y*exp(3*x) + 2*cos(2*y)*cos(3*x) - 4
      u_xy = 6*y*exp(3*x) + 4*y - 6*sin(3*x)*cos(2*y)
      u_xx = 9*y*y*exp(3*x) - 9*sin(2*y)*cos(3*x) + 6
      u_yy = 4*x + 6*y + 2*exp(3*x) - 4*sin(2*y)*cos(3*x)
   end subroutine analytical_test_function_2D


end module FD_test_module

