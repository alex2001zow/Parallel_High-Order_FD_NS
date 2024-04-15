program main

   call run_simulation()

end program main

subroutine run_simulation()
   use omp_lib
   use utility_functions_module, only: read_input_from_command_line
   use mpi_wrapper_module, only: initialize_mpi_wrapper, finalize_mpi_wrapper

   use poisson_module, only: Poission_1D_analytical, Poission_2D_analytical, Poission_3D_analytical
   implicit none

   integer :: ndims, num_physical_cores
   integer, dimension(:), allocatable :: grid_size, processor_dims, stencil_sizes
   real, dimension(:), allocatable :: domain_begin, domain_end

   integer :: rank, world_size

   ! Set the number of threads to the number of physical cores. Hopefully no hyperthreading
   num_physical_cores = omp_get_num_procs()
   call omp_set_num_threads(num_physical_cores)

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size)

   ! MPI-setup. We can also read input from the command line
   ! call read_input_from_command_line(ndims, grid_size, processor_dim)
   ndims = 2

   allocate(grid_size(ndims))
   allocate(processor_dims(ndims))
   allocate(domain_begin(ndims))
   allocate(domain_end(ndims))
   allocate(stencil_sizes(ndims))

   grid_size = 16
   processor_dims = 1
   domain_begin = 0
   domain_end = 1
   stencil_sizes = 3

   if(ndims == 1) then
      call Poission_1D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, stencil_sizes)
   end if
   if(ndims == 2) then
      call Poission_2D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, stencil_sizes)
   end if
   if(ndims == 3) then
      call Poission_3D_analytical(rank, world_size, grid_size, processor_dims, domain_begin, domain_end, stencil_sizes)
   end if

   ! Finalize MPI
   call finalize_mpi_wrapper()

   ! Deallocate the system setup
   if (allocated(grid_size)) deallocate(grid_size)
   if (allocated(processor_dims)) deallocate(processor_dims)
   if (allocated(domain_begin)) deallocate(domain_begin)
   if (allocated(domain_end)) deallocate(domain_end)

end subroutine run_simulation
