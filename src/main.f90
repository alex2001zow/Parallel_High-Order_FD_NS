program main

   call run_simulation()

end program main

subroutine run_simulation()
   use omp_lib
   use utility_functions_module, only: read_input_from_command_line
   use mpi_wrapper_module, only: initialize_mpi_wrapper, finalize_mpi_wrapper

   use FD_test_module, only: FD_test_main
   use block_test_module, only: block_test_main
   use test_Poisson_module, only: Poisson_main
   use Lid_driven_cavity_benchmark_module, only: Lid_driven_cavity_benchmark_2D
   use TravellingWave_2D_module, only: TravelingWave_Poisson_2D_main
   implicit none

   integer :: num_physical_cores, rank, world_size

   ! Set the number of threads to the number of physical cores. Hopefully no hyperthreading.
   num_physical_cores = omp_get_num_procs()
   call omp_set_num_threads(num_physical_cores)

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size)

   ! Test cases
   !call FD_test_main(rank, world_size)
   !call block_test_main(rank, world_size)

   ! Main functions
   !call Poisson_main(rank, world_size)
   call Lid_driven_cavity_benchmark_2D(rank, world_size)
   !call TravelingWave_Poisson_2D_main(rank, world_size)

   ! Finalize MPI
   call finalize_mpi_wrapper()

end subroutine run_simulation
