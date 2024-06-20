module scalapack_module
   use mpi_wrapper_module, only: initialize_mpi_wrapper, finalize_mpi_wrapper
   implicit none

   private

   public :: solve_pde_with_scalapack

contains

   subroutine solve_pde_with_scalapack(iam, numprocs)

      integer, intent(in) :: iam, numprocs
      integer, parameter :: n = 4   ! Size of the matrix
      integer :: nprow, npcol, myrow, mycol, myrank
      integer :: context, info, ictxt, nprocs
      integer :: descA(9), descB(9), ipiv(n)
      double precision :: A(n, n), b(n), x(n)
      double precision, allocatable :: localA(:,:), localB(:)

      ! Initialize the process grid
      nprow = 1
      npcol = 1
      call BLACS_GET(0, 0, context)
      call BLACS_GRIDINIT(context, 'Col', nprow, npcol)
      call BLACS_GRIDINFO(context, nprow, npcol, myrow, mycol)

      ! Define local matrix sizes for each process
      allocate(localA(n/nprow, n/npcol))
      allocate(localB(n/nprow))

      ! Initialize matrix A and vector b on the root process
      if (iam == 0) then
         A = reshape((/ 1.0, 2.0, 0.0, 4.0, &
            5.0, 6.0, 7.0, 8.0, &
            9.0, 10.0, 11.0, 12.0, &
            13.0, 14.0, 15.0, 16.0 /), shape(A))
         b = (/ 1.0, 1.0, 1.0, 1.0 /)
      end if

      ! Descriptor for the distributed matrix A
      call DESCINIT(descA, n, n, n/nprow, n/npcol, 0, 0, context, n/nprow, info)
      if (info /= 0) print *, 'DESCINIT for A returned info = ', info
      call DESCINIT(descB, n, 1, n/nprow, 1, 0, 0, context, n/nprow, info)
      if (info /= 0) print *, 'DESCINIT for B returned info = ', info

      ! Print descriptors
      if (iam == 0) then
         print *, 'Descriptor A:', descA
         print *, 'Descriptor B:', descB
      end if

      ! Distribute matrix A and vector b to the processes
      call pdgemr2d(n, n, A, 1, 1, descA, localA, 1, 1, descA, context)
      call pdgemr2d(n, 1, b, 1, 1, descB, localB, 1, 1, descB, context)

      ! Solve the system Ax = b using ScaLAPACK
      call pdgesv(n, 1, localA, 1, 1, descA, ipiv, localB, 1, 1, descB, info)
      if (info /= 0) print *, 'PDGESV returned info = ', info

      ! Gather the solution vector x back to the root process
      call pdgemr2d(n, 1, localB, 1, 1, descB, b, 1, 1, descB, context)

      ! Output the solution on the root process
      if (iam == 0) then
         print *, 'Solution vector x:'
         print *, b
      end if

      ! Finalize the process grid and MPI
      call blacs_gridexit(context)

   end subroutine solve_pde_with_scalapack

end module scalapack_module

