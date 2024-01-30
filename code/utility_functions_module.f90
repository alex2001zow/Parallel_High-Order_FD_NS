module utility_functions_module
    implicit none
  
    private
    public :: IDX_2D, IDX_3D, print_cartesian_grid
  
    contains

    !> Global index from cartesian system. SHOULD BE COLUMN MAJOR! WE SHOULD ALSO START AT 1. NOT SURE IF CORRECT!
    !! We subtract 1 from ii to start from zero.
    function IDX_2D(N, ii, jj)
      integer, intent(in) :: N, ii, jj
      integer :: IDX_2D
      
      IDX_2D = N * (ii-1) + jj
    end function IDX_2D

    !> Global index from cartesian system. SHOULD BE COLUMN MAJOR! WE SHOULD ALSO START AT 1. NOT SURE IF CORRECT!
    !! We subtract 1 from ii and jj to start from zero.
    function IDX_3D(N, K, ii, jj, kk)
      integer, intent(in) :: N, K, ii, jj, kk
      integer :: IDX_3D
      
      IDX_3D = N * K * (ii-1) + N * (jj-1) + kk
    end function IDX_3D
  
    !> A routine to print the cartesian grid. Just for debugging and understanding
    subroutine print_cartesian_grid(ndim, pn)
      implicit none
      integer, intent(in) :: ndim  ! Number of dimensions
      integer, dimension(ndim), intent(in) :: pn  ! Processors in each dimension
      integer :: i, j, k, idx

      print *, "Cartesian processor grid with dimension:", ndim

      select case(ndim)
      case(1)
          ! 1D Grid
          do i = 1, pn(1)
              write(*, '(I4, 1X)', advance="no") i
          end do
          print *

      case(2)
          ! 2D Grid
          do i = 1, pn(1)
              do j = 1, pn(2)
                  idx = IDX_2D(pn(2), i, j-1)
                  write(*, '(I4, 1X)', advance="no") idx
              end do
              print *  ! New line for the next row
          end do

      case(3)
          ! 3D Grid (printed as slices of 2D grids)
          do i = 1, pn(1)
              print *, "Slice", k, ":"
              do j = 1, pn(2)
                  do k = 1, pn(3)
                      idx = IDX_3D(pn(2), pn(3), i, j, k-1)
                      write(*, '(I4, 1X)', advance="no") idx
                  end do
                  print *  ! New line for the next row
              end do
              if (k < pn(3)) then
                  print *  ! Separate slices with a blank line for readability
              endif
          end do

      end select

  end subroutine print_cartesian_grid
  
end module utility_functions_module
  