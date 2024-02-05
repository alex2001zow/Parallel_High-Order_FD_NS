program main
   use constants_module
   use mpi_wrapper_module
   use rank_parameters_module
   use utility_functions_module
   implicit none

   type(rank_struct) :: rank_params

   integer :: rank, world_size, ierr

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size, ierr)

   ! Setup rank parameters
   call setup_rank_parameters( rank, world_size, rank_params)

   ! Print out rank parameters from rank 0
   if(rank_params%rank == MASTER_RANK) THEN
      call print_cartesian_grid(rank_params%ndims, rank_params%processor_dim)
      call print_rank_parameters(rank_params)
   END IF

   ! Deallocate rank parameters
   call deallocate_rank_parameters(rank_params)

   ! Finalize MPI
   call finalize_mpi_wrapper(ierr)
end program main

