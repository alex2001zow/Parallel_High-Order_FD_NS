module block_test_module
   use constants_module, only : pi, MASTER_RANK
   use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_WTIME
   use mpi_wrapper_module, only: all_reduce_mpi_wrapper, write_block_data_to_file

   use comm_module, only: comm_type, create_cart_comm_type, deallocate_cart_comm_type, &
      print_cart_comm_type, print_cartesian_grid
   use FD_module, only: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil, apply_FDstencil
   use block_module, only: block_type, create_block_type, deallocate_block_type, print_block_type, sendrecv_data_neighbors
   use functions_module, only: FunctionPair, set_function_pointers, calculate_point

   use utility_functions_module, only: open_txt_file, close_txt_file, IDX_XD, sleeper_function
   use initialization_module, only: write_initial_condition_and_boundary

   use solver_module, only: SolverPtrType, set_SystemSolver_pointer, SolverParamsType, set_SolverParamsType, &
      ResultType, print_resultType, choose_iterative_solver
   implicit none

   private

   public :: block_test_main

contains

   subroutine block_test_main(rank, world_size)
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
      processor_dims = 2
      domain_begin = 0
      domain_end = 1
      stencil_sizes = 3

      ! Set the solver parameters tol, max_iter, Jacobi=1 and GS=2
      call set_SolverParamsType(1e-6, 10000, 2, solver_params)

      call block_test(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         stencil_sizes, solver_params)

      ! Deallocate the system setup
      if (allocated(grid_size)) deallocate(grid_size)
      if (allocated(processor_dims)) deallocate(processor_dims)
      if (allocated(domain_begin)) deallocate(domain_begin)
      if (allocated(domain_end)) deallocate(domain_end)

   end subroutine block_test_main

   subroutine block_test(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      stencil_sizes, solver_params)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(:), intent(in) :: grid_size, processor_dims, stencil_sizes
      real, dimension(:), intent(in) :: domain_begin, domain_end
      type(SolverParamsType), intent(in) :: solver_params

      integer, parameter :: num_derivatives = 2
      integer, dimension(num_derivatives), parameter :: derivatives = [1,2]

      call block_test_test(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
         num_derivatives, derivatives, stencil_sizes, solver_params)

   end subroutine block_test

   subroutine block_test_test(ndims, rank, world_size, grid_size, processor_dims, domain_begin, domain_end, &
      num_derivatives, derivatives, stencil_sizes, solver_params)
      integer, intent(in) :: ndims, rank, world_size
      integer, dimension(:), intent(in) :: grid_size, processor_dims
      real, dimension(:), intent(in) :: domain_begin, domain_end

      integer, intent(in) :: num_derivatives
      integer, dimension(:), intent(in) :: derivatives, stencil_sizes
      type(SolverParamsType), intent(in) :: solver_params

      integer, dimension(5) :: num_data_elements

      integer, dimension(ndims) ::  begin, end, temp_ghost_size, temp_stencil_size
      integer :: iounit
      real, dimension(4) :: result_array_with_timings

      type(comm_type) :: comm_params
      type(block_type) :: block_params

      num_data_elements = 1

      call create_cart_comm_type(ndims, processor_dims, rank, world_size, comm_params)

      temp_ghost_size = 2
      temp_stencil_size = stencil_sizes

      !call sleeper_function(1)

      call create_block_type(ndims, num_data_elements(1), num_data_elements(1), temp_ghost_size, temp_stencil_size, &
         domain_begin, domain_end, grid_size, comm_params, block_params)

      block_params%matrix = rank

      !print *, "Before sendrecv: ", "Rank: ", rank, " has the following matrix: ", block_params%matrix

      call sendrecv_data_neighbors(comm_params%comm, block_params)

      !call sendrecv_data_neighbors(comm_params%comm, block_params)

      !print *, "After sendrecv: ", "Rank: ", rank, " has the following matrix: ", block_params%matrix

      ! Write out our setup to a file. Just for debugging purposes
      call open_txt_file("output/output_from_", rank, iounit)

      call print_cart_comm_type(comm_params, iounit)

      call print_block_type(ndims, block_params, iounit)

      call close_txt_file(iounit)

      ! Write out system solution to a file
      call write_block_data_to_file(block_params%data_layout, "output/system_solution.dat", comm_params%comm, block_params%matrix)

      ! Deallocate data

      call deallocate_cart_comm_type(comm_params)

      call deallocate_block_type(block_params)

   end subroutine block_test_test

end module block_test_module





