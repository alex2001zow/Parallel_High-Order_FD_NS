program main
   use mpi_wrapper_module
   use rank_parameters_module
   use utility_functions_module
   implicit none

   type(rank_struct) :: rank_params

   integer :: rank, world_size, ierr

   real, allocatable :: local_matrix(:)

   ! Initialize MPI
   call initialize_mpi_wrapper(rank, world_size, ierr)

   ! Setup rank parameters
   call setup_rank_parameters(rank_params, rank, world_size)

   ! Print out rank parameters from rank 0
   if(rank_params%rank == 0) THEN
      call print_cartesian_grid(rank_params%ndims, rank_params%processor_dim)
      call print_rank_parameters(rank_params)
   END IF

   allocate(local_matrix(rank_params%total_num_elements))
   deallocate(local_matrix)

   ! Deallocate rank parameters
   call deallocate_rank_parameters(rank_params)

   ! Finalize MPI
   call finalize_mpi_wrapper(ierr)
end program main

