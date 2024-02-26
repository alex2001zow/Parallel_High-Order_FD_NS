module finite_difference_module
   use utility_functions_module, only : IDX_XD_REVERSED, sleeper_function, print_matrix, flip_1d_integer_array
   use block_module, only : block_type
   implicit none

   private

   type FDstencil_type
      integer :: num_derivatives
      integer, dimension(:), allocatable :: derivatives, alphas, betas, stencil_sizes
      real, dimension(:), allocatable :: derivatives_sign, dx, stencil_coefficients
   end type FDstencil_type

   public :: FDstencil_type, create_finite_difference_stencil, deallocate_finite_difference_stencil,&
      print_finite_difference_stencil, apply_FDstencil

contains

   !> Create a finite difference stencil
   subroutine create_finite_difference_stencil(ndims, num_derivatives, derivatives, derivatives_sign, &
      dx, alphas, betas, FDstencil_type_output)
      integer, intent(in) :: ndims, num_derivatives
      integer, dimension(ndims * num_derivatives), intent(in) :: derivatives
      integer, dimension(ndims), intent(in) :: alphas, betas
      real,  dimension(ndims), intent(in) :: dx
      real,  dimension(num_derivatives), intent(in) :: derivatives_sign

      real, dimension(:), allocatable :: temp_stencil_coefficients
      type(FDstencil_type), intent(out) :: FDstencil_type_output

      integer :: ii
      integer, dimension(ndims) :: derivative

      allocate(FDstencil_type_output%derivatives(ndims * num_derivatives))
      allocate(FDstencil_type_output%derivatives_sign(num_derivatives))
      allocate(FDstencil_type_output%dx(ndims))
      allocate(FDstencil_type_output%alphas(ndims))
      allocate(FDstencil_type_output%betas(ndims))
      allocate(FDstencil_type_output%stencil_sizes(ndims))

      FDstencil_type_output%num_derivatives = num_derivatives
      FDstencil_type_output%derivatives = derivatives
      FDstencil_type_output%derivatives_sign = derivatives_sign
      FDstencil_type_output%dx = dx
      FDstencil_type_output%alphas = alphas
      FDstencil_type_output%betas = betas
      FDstencil_type_output%stencil_sizes = alphas + betas + 1

      allocate(FDstencil_type_output%stencil_coefficients(product(FDstencil_type_output%stencil_sizes)))
      FDstencil_type_output%stencil_coefficients(:) = 0.0

      allocate(temp_stencil_coefficients(product(FDstencil_type_output%stencil_sizes)))

      do ii = 1, num_derivatives
         derivative = derivatives((ii-1)*ndims+1:ii*ndims)
         call finite_difference_stencil(ndims, FDstencil_type_output%alphas, FDstencil_type_output%betas, &
            FDstencil_type_output%stencil_sizes, derivative, FDstencil_type_output%dx, &
            temp_stencil_coefficients)

         FDstencil_type_output%stencil_coefficients = FDstencil_type_output%stencil_coefficients + &
            derivatives_sign(ii) * temp_stencil_coefficients
      end do

      deallocate(temp_stencil_coefficients)

   end subroutine create_finite_difference_stencil

   !> Destroy a finite difference stencil
   subroutine deallocate_finite_difference_stencil(FDstencil_type_input)
      type(FDstencil_type), intent(inout) :: FDstencil_type_input

      if(allocated(FDstencil_type_input%derivatives)) deallocate(FDstencil_type_input%derivatives)
      if(allocated(FDstencil_type_input%derivatives_sign)) deallocate(FDstencil_type_input%derivatives_sign)
      if(allocated(FDstencil_type_input%dx)) deallocate(FDstencil_type_input%dx)
      if(allocated(FDstencil_type_input%alphas)) deallocate(FDstencil_type_input%alphas)
      if(allocated(FDstencil_type_input%betas)) deallocate(FDstencil_type_input%betas)
      if(allocated(FDstencil_type_input%stencil_sizes)) deallocate(FDstencil_type_input%stencil_sizes)

      if(allocated(FDstencil_type_input%stencil_coefficients)) deallocate(FDstencil_type_input%stencil_coefficients)

   end subroutine deallocate_finite_difference_stencil

   !> Print the finite difference stencil
   subroutine print_finite_difference_stencil(ndims, FDstencil_type_input, iounit)
      integer, intent(in) :: ndims, iounit
      type(FDstencil_type), intent(in) :: FDstencil_type_input

      integer :: ii, jj, global_index
      integer, dimension(ndims) :: index

      index(:) = 1

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "FDstencil_type: "
      write(iounit, *)

      write(iounit, *) "num_derivatives: ", FDstencil_type_input%num_derivatives
      write(iounit, *) "derivatives: ", FDstencil_type_input%derivatives
      write(iounit, '(A, *(F10.3))') "derivatives_sign: ", FDstencil_type_input%derivatives_sign
      write(iounit, '(A, *(F10.3))') "dx: ", FDstencil_type_input%dx
      write(iounit, *) "alphas: ", FDstencil_type_input%alphas
      write(iounit, *) "betas: ", FDstencil_type_input%betas
      write(iounit, *) "stencil_sizes: ", FDstencil_type_input%stencil_sizes

      write(iounit, *)
      write(iounit, *) "stencil_coefficients: "
      if(ndims == 2) then
         do ii = -FDstencil_type_input%alphas(1), FDstencil_type_input%betas(1)
            do jj = -FDstencil_type_input%alphas(2), FDstencil_type_input%betas(2)
               global_index = IDX_XD_REVERSED(ndims, FDstencil_type_input%stencil_sizes, index)
               if (ii == 0 .and. jj == 0) then
                  write(iounit, '(F10.3,A)', advance='no') FDstencil_type_input%stencil_coefficients(global_index), "*"
               else
                  write(iounit, '(F10.3)', advance='no') FDstencil_type_input%stencil_coefficients(global_index)
               endif
               index(2) = index(2) + 1
            end do
            write(iounit, *)
            index(1) = index(1) + 1
            index(2) = 1
         end do
      end if

      if ( ndims == 3 ) then
         write(iounit, *) "NOT DONE YET"
      end if

   end subroutine print_finite_difference_stencil

   !> Calculate the finite difference stencil coefficients.
   !! sum_{n=-alpha_1}^{beta_1} sum_{m=-alpha_2}^{beta_2} a_nm * A_ij
   subroutine finite_difference_stencil(ndims, alphas, betas, stencil_sizes, derivative, dx, stencil_coefficients)
      integer, intent(in) :: ndims

      integer, dimension(ndims), intent(in) :: alphas, betas, stencil_sizes, derivative
      real, dimension(ndims), intent(in) :: dx

      real, dimension(product(stencil_sizes)), intent(out) :: stencil_coefficients

      integer :: ii, jj, global_index, stride, info
      integer, dimension(product(stencil_sizes)) :: ipiv  ! 'ipiv' is an array for pivot indices used by DGESV
      real, dimension(product(stencil_sizes)) :: taylor_coefficients
      real, dimension((product(stencil_sizes)**ndims)) :: coefficient_matrix

      info = 0

      stride = product(stencil_sizes)

      ! Initialize the stencil coefficients. It starts as the rhs of the linear system.
      stencil_coefficients = 0.0
      stencil_coefficients(IDX_XD_REVERSED(ndims, stencil_sizes, derivative)) = 1.0

      global_index = 1
      do ii = -alphas(1), betas(1)
         do jj = -alphas(2), betas(2)
            call calculate_taylor_coefficients(ndims, stencil_sizes, dx * [ii, jj], taylor_coefficients)

            coefficient_matrix(global_index:global_index + stride - 1) = taylor_coefficients
            global_index = global_index + stride
         end do
      end do

      !Call DGESV to solve the linear system
      call DGESV(stride, 1, coefficient_matrix, stride, ipiv, stencil_coefficients, stride, info)

      if (info /= 0) then
         print *, "DGESV reported an error: ", info
         stop
      end if

   end subroutine finite_difference_stencil

   !> Calculate the Taylor coefficients up to order s
   !! A_ij = sum_{i=0}^{p-1} sum_{j=0}^{q-1} 1/(i! * j!) (n*h)^i * (m*k)^j
   !! Here the factorial is slow. We could use the loop to calculate the factorial. But this is for another time.
   subroutine calculate_taylor_coefficients(ndims, stencil_dims, dx, taylor_coefficients)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: stencil_dims
      real, dimension(ndims), intent(in) :: dx
      real, dimension(product(stencil_dims)), intent(out) :: taylor_coefficients

      integer :: ii, jj, global_index
      real :: numerator(ndims), denominator(ndims)

      numerator = 0.0
      denominator = 0.0
      taylor_coefficients = 0.0

      do ii = 0, stencil_dims(1)-1
         numerator(1) = dx(1) ** ii
         denominator(1) = gamma(real(ii+1))
         do jj = 0, stencil_dims(2)-1
            numerator(2) = dx(2) ** jj
            denominator(2) = gamma(real(jj+1))
            global_index = IDX_XD_REVERSED(ndims, stencil_dims, [ii+1,jj+1])
            taylor_coefficients(global_index) = product(numerator) / product(denominator)
         end do
      end do

   end subroutine calculate_taylor_coefficients

   !> Apply the finite difference stencil to a matrix
   !! Replace IDX_XD with a more efficient way to calculate the global index
   function apply_FDstencil(ndims, FDstencil_type_input, block, index) result(value)
      integer, intent(in) :: ndims
      type(FDstencil_type), intent(in) :: FDstencil_type_input
      type(block_type), intent(in) :: block
      integer, dimension(ndims), intent(in) :: index

      real :: value
      integer :: ii, jj, global_index_stencil, global_index_matrix

      integer, dimension(ndims) :: index_stencil

      index_stencil = 1

      value = 0.0

      do ii = -FDstencil_type_input%alphas(1), FDstencil_type_input%betas(1)
         do jj = -FDstencil_type_input%alphas(2), FDstencil_type_input%betas(2)
            global_index_matrix = IDX_XD_REVERSED(ndims, block%size, index + [ii,jj])
            global_index_stencil = IDX_XD_REVERSED(ndims, FDstencil_type_input%stencil_sizes, index_stencil)
            value = value + FDstencil_type_input%stencil_coefficients(global_index_stencil) * block%matrix(global_index_matrix)
            index_stencil(2) = index_stencil(2) + 1
         end do
         index_stencil(1) = index_stencil(1) + 1
         index_stencil(2) = 1
      end do

   end function apply_FDstencil

end module finite_difference_module
