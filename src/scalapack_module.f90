module scalapack_module
   implicit none
   private
   public :: solve_pde_with_scalapack, solve_pde_with_scalapack_2

contains
   subroutine solve_pde_with_scalapack(rank, world_size)
      integer, intent(in) :: rank, world_size
      integer, parameter :: n = 4 ! Size of the matrix
      integer :: nprow, npcol, myrow, mycol, myrank
      integer :: context, info, ictxt, nprocs
      integer :: descA(9), descB(9)
      double precision :: A(n, n), b(n), x(n), work(n)
      integer :: ipiv(n)

      ! Initialize the process grid
      nprow = 1
      npcol = 1
      call BLACS_GET(0, 0, context)
      call BLACS_GRIDINIT(context, 'R', nprow, npcol)
      call BLACS_GRIDINFO(context, nprow, npcol, myrow, mycol)

      ! Initialize matrix A and vector b on the root process
      if (rank == 0) then
         A = reshape((/ 1.0, 0.0, 0.0, 0.0, &
            0.0, 2.0, 0.0, 0.0, &
            0.0, 0.0, 3.0, 0.0, &
            0.0, 0.0, 0.0, 4.0 /), shape(A))
         b = (/ 1.0, 2.0, 3.0, 4.0 /)
      end if

      ! Descriptor for the distributed matrix A
      call DESCINIT(descA, n, n, n, n, 0, 0, context, n, info)
      if (info /= 0) print *, 'DESCINIT for A returned info = ', info

      call DESCINIT(descB, n, 1, n, 1, 0, 0, context, n, info)
      if (info /= 0) print *, 'DESCINIT for B returned info = ', info

      ! Solve the system Ax = b using ScaLAPACK
      call pdgesv(n, 1, A, 1, 1, descA, ipiv, b, 1, 1, descB, info)
      if (info /= 0) print *, 'PDGESV returned info = ', info

      ! Output the solution on the root process
      if (rank == 0) then
         print *, 'Solution vector x:'
         print *, b
      end if

      ! Finalize the process grid
      call blacs_gridexit(context)
   end subroutine solve_pde_with_scalapack

   subroutine solve_pde_with_scalapack_2(rank, size)
      integer, intent(in) :: rank, size

      integer, parameter :: n = 4 ! Size of the matrix
      integer :: nprow, npcol, myrow, mycol, info
      integer :: context, ictxt
      integer :: descA(9), descB(9)
      double precision, allocatable :: A(:,:), b(:)
      integer, allocatable :: ipiv(:)
      integer :: nb = 64

      if (size /= 2) then
         if (rank == 0) print *, "This program should be run with exactly 2 processes."
         return
      end if

      ! Set up process grid (1x2)
      nprow = 1
      npcol = 2

      ! Initialize BLACS
      call BLACS_GET(-1, 0, context)
      call BLACS_GRIDINIT(context, 'C', nprow, npcol)
      call BLACS_GRIDINFO(context, nprow, npcol, myrow, mycol)

      ! Allocate local arrays
      allocate(A(n/npcol, n/nprow))
      allocate(b(n/npcol))
      allocate(ipiv(n/npcol + nb))  ! nb is the blocking factor, often 64

      ! Initialize matrix A and vector b on all processes
      if (mycol == 0) then
         A = reshape((/ 1.0, 2.0, 5.0, 6.0, 9.0, 10.0, 13.0, 14.0/), shape(A))
         b = (/ 1.0, 2.0 /)
      else
         A = reshape((/ 3.0, 4.0, 7.0, 8.0, 11.0, 12.0, 15.0, 16.0/), shape(A))
         b = (/ 3.0, 4.0 /)
      end if

      ! Set up array descriptors
      call DESCINIT(descA, n, n, n/npcol, n/nprow, 0, 0, context, n, info)
      if (info /= 0) print *, 'DESCINIT for A returned info = ', info

      call DESCINIT(descB, n, 1, n/npcol, 1, 0, 0, context, n, info)
      if (info /= 0) print *, 'DESCINIT for B returned info = ', info

      ! Solve the system Ax = b using ScaLAPACK
      call PDGESV(n, 1, A, 1, 1, descA, ipiv, b, 1, 1, descB, info)
      if (info /= 0) print *, 'PDGESV returned info = ', info

      ! Gather the solution to rank 0 and print
      if (rank == 0) then
         call PDGEMR2D(n, 1, b, 1, 1, descB, b, 1, 1, descB, context)
         print *, 'Solution vector x:'
         print *, b
      end if

      ! Clean up
      deallocate(A, b, ipiv)
      call BLACS_GRIDEXIT(context)

   end subroutine solve_pde_with_scalapack_2

end module scalapack_module

! 1 2 3 4
! 5 6 7 8
! 9 10 11 12
! 13 14 15 16

! 1
! 2
! 3
! 4
