module finite_difference_module
   use utility_functions_module, only : IDX_XD
   implicit none

   private

   type FDstencil_type
      integer, dimension(:), allocatable ::derivatives, alphas, betas, stencil_sizes
      real, dimension(:), allocatable :: dx
      real, dimension(:), allocatable :: stencil_coefficients
   end type FDstencil_type

   public :: FDstencil_type, create_finite_difference_stencil, deallocate_finite_difference_stencil,&
      print_finite_difference_stencil

contains

   !> Create a finite difference stencil
   subroutine create_finite_difference_stencil(ndims, derivatives, dx, alphas, betas, FDstencil_type_output)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: derivatives, alphas, betas
      real,  dimension(ndims), intent(in) :: dx

      type(FDstencil_type), intent(out) :: FDstencil_type_output

      allocate(FDstencil_type_output%derivatives(ndims))
      allocate(FDstencil_type_output%alphas(ndims))
      allocate(FDstencil_type_output%betas(ndims))
      allocate(FDstencil_type_output%stencil_sizes(ndims))
      allocate(FDstencil_type_output%dx(ndims))

      FDstencil_type_output%derivatives = derivatives
      FDstencil_type_output%alphas = alphas
      FDstencil_type_output%betas = betas
      FDstencil_type_output%stencil_sizes = alphas + betas + 1
      FDstencil_type_output%dx = dx

      allocate(FDstencil_type_output%stencil_coefficients(product(FDstencil_type_output%stencil_sizes)))

      call finite_difference_stencil(ndims, FDstencil_type_output)

   end subroutine create_finite_difference_stencil

   !> Destroy a finite difference stencil
   subroutine deallocate_finite_difference_stencil(FDstencil_type_input)
      type(FDstencil_type), intent(inout) :: FDstencil_type_input

      deallocate(FDstencil_type_input%derivatives)
      deallocate(FDstencil_type_input%alphas)
      deallocate(FDstencil_type_input%betas)
      deallocate(FDstencil_type_input%stencil_sizes)
      deallocate(FDstencil_type_input%dx)

      deallocate(FDstencil_type_input%stencil_coefficients)

   end subroutine deallocate_finite_difference_stencil

   !> Print the finite difference stencil
   !! stencil: the stencil to print
   subroutine print_finite_difference_stencil(ndims, FDstencil_type_input, iounit)
      integer, intent(in) :: ndims, iounit
      type(FDstencil_type), intent(in) :: FDstencil_type_input

      integer :: ii, jj, global_index

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "FDstencil_type: "
      write(iounit, *)

      write(iounit, *) "derivatives: ", FDstencil_type_input%derivatives
      write(iounit, *) "alphas: ", FDstencil_type_input%alphas
      write(iounit, *) "betas: ", FDstencil_type_input%betas
      write(iounit, *) "stencil_sizes: ", FDstencil_type_input%stencil_sizes
      write(iounit, '(A, *(F10.3))') "dx: ", FDstencil_type_input%dx

      write(iounit, *)
      write(iounit, *) "stencil_coefficients: "
      if(ndims == 2) then

         do ii = 1, FDstencil_type_input%stencil_sizes(1)
            do jj = 1, FDstencil_type_input%stencil_sizes(2)
               global_index = IDX_XD(2, FDstencil_type_input%stencil_sizes, [jj,ii])
               write(iounit, '(F10.3)', advance='no') FDstencil_type_input%stencil_coefficients(global_index)
            end do
            write(iounit, *)
         end do

      end if

      if ( ndims == 3 ) then
         write(iounit, *) "NOT DONE YET"
      end if

   end subroutine print_finite_difference_stencil

   !> Calculate the finite difference stencil coefficients.
   !! sum_{n=-alpha}^{beta} a_n * sum_{k=0}^{s = alpha + beta + 1} d^k/dx^k f(x)/k! *(n*h)^k
   subroutine finite_difference_stencil(ndims, FDstencil_type_inout)
      integer, intent(in) :: ndims
      type(FDstencil_type), intent(inout) :: FDstencil_type_inout

      integer :: ii, jj, global_index, stride, info
      integer, dimension(ndims) :: index
      integer, dimension(product(FDstencil_type_inout%stencil_sizes)) :: ipiv  ! 'ipiv' is an array for pivot indices used by DGESV
      real, dimension(product(FDstencil_type_inout%stencil_sizes)) :: taylor_coefficients
      real, dimension((product(FDstencil_type_inout%stencil_sizes))**ndims) :: coefficient_matrix

      info = 0

      stride = product(FDstencil_type_inout%stencil_sizes)

      ! Initialize the stencil coefficients. It starts as the rhs of the linear system.
      FDstencil_type_inout%stencil_coefficients(:) = 0.0
      FDstencil_type_inout%stencil_coefficients(IDX_XD(ndims,FDstencil_type_inout%stencil_sizes,&
         FDstencil_type_inout%derivatives)) = 1.0

      index = [1,1]
      do ii = -FDstencil_type_inout%alphas(1), FDstencil_type_inout%betas(1)
         index(1) = 1
         do jj = -FDstencil_type_inout%alphas(2), FDstencil_type_inout%betas(2)
            ! Calculate the global index for the coefficient matrix
            global_index = IDX_XD(2, FDstencil_type_inout%stencil_sizes, index)
            call calculate_taylor_coefficients(ndims, FDstencil_type_inout%stencil_sizes, &
               FDstencil_type_inout%dx, index, taylor_coefficients)

            coefficient_matrix(global_index:global_index + stride - 1) = taylor_coefficients
            index(1) = index(1) + 1
         end do
         index(2) = index(2) + 1
      end do

      ! Call DGESV to solve the linear system
      !call DGESV(stride, 1, coefficient_matrix, stride, ipiv, FDstencil_type_inout%stencil_coefficients, stride, info)

      if (info /= 0) then
         print *, "DGESV reported an error: ", info
         stop
      end if

   end subroutine finite_difference_stencil

   !> Calculate the Taylor coefficients up to order s
   !! sum_{i=0}^{p-1} sum_{j=0}^{q-1} 1/(i! * j!) (n*h)^i * (m*k)^j
   !! Double check that ii and jj are correct. It should be running from 0 to p-1 and 0 to q-1
   subroutine calculate_taylor_coefficients(ndims, stencil_sizes, dx, index, coefficients)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: stencil_sizes, index
      real, dimension(ndims), intent(in) :: dx
      real, dimension(product(stencil_sizes)), intent(out) :: coefficients

      integer :: ii, jj, global_index
      real, dimension(ndims) :: numerator, denominator

      numerator(:) = 1.0
      denominator(:) = 1.0

      if(ndims == 2) then
         do ii = 1, stencil_sizes(1)
            numerator(2) = 1.0
            denominator(2) = 1.0
            do jj = 1, stencil_sizes(2)
               global_index = IDX_XD(2, stencil_sizes, [jj,ii])
               coefficients(global_index) = product(numerator) / product(denominator)
               numerator(2) = numerator(2) * index(2) * dx(2)
               denominator(2) = denominator(2) * jj
            end do
            numerator(1) = numerator(1) * index(1) * dx(1)
            denominator(1) = denominator(1) * ii
         end do
      end if

      if(ndims == 3) then
         print *, "3D not implemented yet"
         stop
      end if

   end subroutine calculate_taylor_coefficients

end module finite_difference_module
