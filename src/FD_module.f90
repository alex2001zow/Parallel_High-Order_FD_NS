module FD_module
   use utility_functions_module, only : IDX_XD, IDX_XD_INV, print_real_array, sleeper_function
   use block_module, only : block_type
   implicit none

   private

   type FDstencil_type
      integer :: num_derivatives
      integer, dimension(:), allocatable :: derivatives, alphas, betas, stencil_sizes
      real, dimension(:), allocatable :: derivatives_sign, dx, stencil_coefficients
      real :: center_coefficient
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

      integer :: ii, index_start, index_end
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
         index_start = (ii-1)*ndims+1
         index_end = ii*ndims
         derivative = derivatives(index_start:index_end)
         call calculate_finite_difference_stencil(ndims, FDstencil_type_output%alphas, FDstencil_type_output%betas, &
            FDstencil_type_output%stencil_sizes, derivative, temp_stencil_coefficients)

         FDstencil_type_output%stencil_coefficients = FDstencil_type_output%stencil_coefficients + &
            derivatives_sign(ii) * temp_stencil_coefficients / product(dx**derivatives(index_start:index_end))
      end do

      !> The center coefficient is the coefficient of the central point of the stencil. Found at alphas + 1. Used to "normalize" the stencil.
      FDstencil_type_output%center_coefficient = FDstencil_type_output%stencil_coefficients(IDX_XD(ndims, &
         FDstencil_type_output%stencil_sizes, alphas + 1))

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

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "FDstencil_type: "
      write(iounit, *)

      write(iounit, *) "num_derivatives: ", FDstencil_type_input%num_derivatives
      write(iounit, *) "derivatives: ", FDstencil_type_input%derivatives
      write(iounit, *) "derivatives_sign: ", FDstencil_type_input%derivatives_sign
      write(iounit, *) "dx: ", FDstencil_type_input%dx
      write(iounit, *) "alphas: ", FDstencil_type_input%alphas
      write(iounit, *) "betas: ", FDstencil_type_input%betas
      write(iounit, *) "stencil_sizes: ", FDstencil_type_input%stencil_sizes
      write(iounit, *) "center_coefficient: ", FDstencil_type_input%center_coefficient

      call print_real_array(ndims, FDstencil_type_input%stencil_sizes, FDstencil_type_input%stencil_coefficients, 1, &
         "stencil_coefficients: ", iounit)

   end subroutine print_finite_difference_stencil

   !> Calculate the finite difference stencil coefficients.
   !! sum_{n=-alpha_1}^{beta_1} sum_{m=-alpha_2}^{beta_2} a_nm * A_ij
   subroutine calculate_finite_difference_stencil(ndims, alphas, betas, stencil_sizes, derivative, stencil_coefficients)
      integer, intent(in) :: ndims

      integer, dimension(ndims), intent(in) :: alphas, betas, stencil_sizes, derivative

      real, dimension(product(stencil_sizes)), intent(out) :: stencil_coefficients

      integer, dimension(ndims) :: index
      integer :: global_index, start_index, end_index, stride, info
      integer, dimension(product(stencil_sizes)) :: ipiv  ! 'ipiv' is an array for pivot indices used by DGESV
      real, dimension(product(stencil_sizes)) :: taylor_coefficients
      real, dimension((product(stencil_sizes*stencil_sizes)**ndims)) :: coefficient_matrix

      info = 0

      stride = product(stencil_sizes)

      ! Initialize the stencil coefficients. It starts as the rhs of the linear system.
      stencil_coefficients = 0.0
      stencil_coefficients(IDX_XD(ndims, stencil_sizes, derivative + 1)) = 1.0 ! derivative + 1 because the index should start at 1.

      do global_index = 1, product(stencil_sizes)
         index = IDX_XD_INV(ndims, stencil_sizes, global_index)
         index = index - alphas - 1 ! go from alpha to beta in steps of 1
         start_index = (global_index - 1) * stride + 1
         end_index = global_index * stride

         ! The Taylor coefficients are calculated for each point of the stencil
         call calculate_taylor_coefficients(ndims, stencil_sizes, real(index,kind=8), taylor_coefficients)
         coefficient_matrix(start_index:end_index) = taylor_coefficients

      end do

      ! Call DGESV to solve the linear system
      call DGESV(stride, 1, coefficient_matrix, stride, ipiv, stencil_coefficients, stride, info)

      if (info /= 0) then
         print *, "DGESV reported an error: ", info
         stop
      end if

   end subroutine calculate_finite_difference_stencil

   !> Calculate the Taylor coefficients up to order s
   !! A_ij = sum_{i=0}^{p-1} sum_{j=0}^{q-1} 1/(i! * j!) (n*h)^i * (m*k)^j
   !! Here the factorial is slow. We could use the loop to calculate the factorial. But this is for another time. Should be fine for setup.
   subroutine calculate_taylor_coefficients(ndims, stencil_sizes, dx, taylor_coefficients)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: stencil_sizes
      real, dimension(ndims), intent(in) :: dx
      real, dimension(product(stencil_sizes)), intent(out) :: taylor_coefficients

      integer :: global_index
      integer, dimension(ndims) :: current_index
      real :: numerator(ndims), denominator(ndims)

      numerator = 0.0
      denominator = 0.0
      taylor_coefficients = 0.0

      do global_index = 1, product(stencil_sizes)
         current_index = IDX_XD_INV(ndims, stencil_sizes, global_index)
         current_index = current_index

         numerator = dx ** real(current_index-1,kind=8)
         denominator = gamma(real(current_index,kind=8))
         taylor_coefficients(global_index) = product(numerator) / product(denominator)
      end do

   end subroutine calculate_taylor_coefficients

   !> Apply the finite difference stencil to a matrix
   !! Replace IDX_XD with a more efficient way to calculate the global index
   function apply_FDstencil(ndims, FDstencil_type_input, dims, matrix, index) result(val)
      integer, intent(in) :: ndims
      type(FDstencil_type), intent(in) :: FDstencil_type_input
      integer, dimension(ndims), intent(in) :: dims
      real, dimension(product(dims)) :: matrix
      integer, dimension(ndims), intent(in) :: index

      real :: stencil_val, matrix_val, val
      integer :: global_index, block_index

      integer, dimension(ndims) :: index_stencil

      val = 0.0

      do global_index = 1, product(FDstencil_type_input%stencil_sizes)
         index_stencil = IDX_XD_INV(ndims, FDstencil_type_input%stencil_sizes, global_index)
         index_stencil = index_stencil - FDstencil_type_input%alphas - 1
         block_index = IDX_XD(ndims, dims, index + index_stencil)

         stencil_val = FDstencil_type_input%stencil_coefficients(global_index)
         matrix_val = matrix(block_index)

         val = val + stencil_val * matrix_val
      end do

      ! Subtract the center coefficient
      val = val - FDstencil_type_input%center_coefficient * matrix(IDX_XD(ndims, dims, index))

   end function apply_FDstencil

end module FD_module
