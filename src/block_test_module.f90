module block_test_module
   use constants_module, only : pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, apply_FDstencil
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type, sendrecv_data_neighbors
   use functions_module, only: FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, IDX_XD_INV, sleeper_function
   use initialization_module, only: write_initial_condition_and_boundary

   use solver_module, only: SolverPtrType, set_SystemSolver_pointer, SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, choose_iterative_solver
   implicit none

   private
   !> Number of dimensions and number of derivatives
   integer, parameter :: ndims = 2

   !> Grid parameters
   integer, dimension(ndims), parameter :: grid_size = 16, processor_dims = [2,1]
   logical, dimension(ndims), parameter :: periods = [.false.,.false.]
   logical, parameter :: reorder = .true.
   real, dimension(ndims), parameter :: domain_begin = 0, domain_end = 1
   integer, dimension(ndims), parameter :: stencil_sizes = 3
   integer, dimension(ndims), parameter :: ghost_begin = [0,0], ghost_end = [0,0]
   integer, dimension(ndims), parameter :: stencil_begin = stencil_sizes/2, stencil_end = stencil_sizes/2

   public :: block_test_main

contains

   subroutine block_test_main(rank, world_size)
      integer, intent(in) :: rank, world_size

      type(comm_type) :: comm_params
      type(block_type) :: block_params

      integer :: iounit

      call create_cart_comm_type(ndims, processor_dims, periods, reorder, rank, world_size, comm_params)

      !call sleeper_function(1)

      call create_block_type(ndims, 1, 1, domain_begin, domain_end, grid_size, comm_params, &
         ghost_begin, ghost_end, stencil_begin, stencil_end, 1, block_params)

      block_params%matrix_ptr = rank

      if(ndims == 1) then
         !call write_global_index_to_block_1D(block_params)
      else if(ndims == 2) then
         !call write_global_index_to_block_2D(block_params)
      end if

      call sendrecv_data_neighbors(comm_params%comm, block_params, block_params%matrix_ptr)

      ! Write out our setup to a file. Just for debugging purposes
      call open_txt_file("output/output_from_", rank, iounit)

      !call print_cart_comm_type(comm_params, iounit)

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

   end subroutine block_test_main

   pure subroutine write_global_index_to_block_1D(block_params)

      type(block_type), intent(inout) :: block_params

      integer :: ii
      integer :: local_index
      integer, dimension(block_params%ndims) :: local_indices

      do ii = block_params%block_begin_c(1)+1, block_params%block_end_c(1)
         local_indices = [ii]
         call IDX_XD(block_params%ndims, block_params%extended_block_dims, local_indices, local_index)
         block_params%matrix_ptr(local_index) = local_index

      end do

   end subroutine write_global_index_to_block_1D

   pure subroutine write_global_index_to_block_2D(block_params)
      type(block_type), intent(inout) :: block_params

      integer :: ii, jj
      integer :: local_index
      integer, dimension(block_params%ndims) :: local_indices

      do ii = block_params%block_begin_c(1)+1, block_params%block_end_c(1)
         do jj = block_params%block_begin_c(2)+1, block_params%block_end_c(2)
            local_indices = [ii, jj]
            call IDX_XD(block_params%ndims, block_params%extended_block_dims, local_indices, local_index)
            block_params%matrix_ptr(local_index) = local_index

         end do
      end do

   end subroutine write_global_index_to_block_2D

end module block_test_module





