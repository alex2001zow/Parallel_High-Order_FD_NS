module FD_module
   use utility_functions_module, only : IDX_XD, IDX_XD_INV, print_real_array, sleeper_function
   use block_module, only : block_type
   implicit none

   private

   type FDstencil_type
      integer :: num_derivatives
      integer, dimension(:), allocatable :: derivatives, derivatives_order, alphas, betas, stencil_sizes
      real, dimension(:), allocatable :: dx, stencil_coefficients
      real :: center_coefficient
   end type FDstencil_type

   public :: FDstencil_type, create_finite_difference_stencil, deallocate_finite_difference_stencil,&
      print_finite_difference_stencil, apply_FDstencil, create_finite_difference_stencil_from_order

contains

   !> subroutine to create a finite difference stencil from a given order. We also try to make a centered stencil first!
   subroutine create_finite_difference_stencil_from_order(ndims, num_derivatives, derivatives, dx, order, &
      begin, end, indices, FDstencil_type_output)
      integer, intent(in) :: ndims, num_derivatives, order
      integer, dimension(ndims * num_derivatives), intent(in) :: derivatives
      real,  dimension(ndims), intent(in) :: dx
      integer, dimension(ndims), intent(in) :: begin, end, indices

      type(FDstencil_type), intent(out) :: FDstencil_type_output

      integer, dimension(ndims) :: dist_to_begin, dist_to_end
      integer, dimension(ndims) :: stencil_sizes, alphas, betas, adjusted_alphas, adjusted_betas, alpha_diff, beta_diff

      ! Find the stencil size that kinda matches the order
      stencil_sizes = find_stencil_size_from_order(ndims, num_derivatives, derivatives, order)

      ! Try with a centered stencil first
      alphas = (stencil_sizes + 1) / 2
      betas = stencil_sizes - alphas

      dist_to_begin = indices - begin
      dist_to_end = end - indices

      ! Adjust alphas and betas to stay within bounds
      adjusted_alphas = max(min(alphas, dist_to_begin), 0)
      adjusted_betas = max(min(betas, dist_to_end), 0)

      ! If one side exceeds the bounds, add the excess to the other side
      alpha_diff = alphas - adjusted_alphas
      beta_diff = betas - adjusted_betas
      adjusted_alphas = adjusted_alphas + beta_diff
      adjusted_betas = adjusted_betas + alpha_diff

      call create_finite_difference_stencil(ndims, num_derivatives, derivatives, dx, &
         adjusted_alphas, adjusted_betas, FDstencil_type_output)

   end subroutine create_finite_difference_stencil_from_order

   !> Create a finite difference stencil
   subroutine create_finite_difference_stencil(ndims, num_derivatives, derivatives, &
      dx, alphas, betas, FDstencil_type_output)
      integer, intent(in) :: ndims, num_derivatives
      integer, dimension(ndims * num_derivatives), intent(in) :: derivatives
      integer, dimension(ndims), intent(in) :: alphas, betas
      real,  dimension(ndims), intent(in) :: dx

      real, dimension(:), allocatable :: temp_stencil_coefficients
      type(FDstencil_type), intent(out) :: FDstencil_type_output

      integer :: global_index, index_start, index_end
      integer, dimension(ndims) :: derivative

      allocate(FDstencil_type_output%derivatives(ndims * num_derivatives))
      allocate(FDstencil_type_output%derivatives_order(num_derivatives))
      allocate(FDstencil_type_output%dx(ndims))
      allocate(FDstencil_type_output%alphas(ndims))
      allocate(FDstencil_type_output%betas(ndims))
      allocate(FDstencil_type_output%stencil_sizes(ndims))

      FDstencil_type_output%num_derivatives = num_derivatives
      FDstencil_type_output%derivatives = derivatives
      FDstencil_type_output%dx = dx
      FDstencil_type_output%alphas = alphas
      FDstencil_type_output%betas = betas
      FDstencil_type_output%stencil_sizes = alphas + betas + 1
      FDstencil_type_output%derivatives_order = calculate_derivative_order(ndims, num_derivatives, derivatives, &
         FDstencil_type_output%stencil_sizes)

      allocate(FDstencil_type_output%stencil_coefficients(product(FDstencil_type_output%stencil_sizes)))
      FDstencil_type_output%stencil_coefficients(:) = 0.0

      allocate(temp_stencil_coefficients(product(FDstencil_type_output%stencil_sizes)))

      do global_index = 1, num_derivatives
         index_start = (global_index-1)*ndims+1
         index_end = global_index*ndims
         derivative = derivatives(index_start:index_end)
         call calculate_finite_difference_stencil(ndims, FDstencil_type_output%alphas, FDstencil_type_output%betas, &
            FDstencil_type_output%stencil_sizes, derivative, temp_stencil_coefficients)

         FDstencil_type_output%stencil_coefficients = FDstencil_type_output%stencil_coefficients + &
            temp_stencil_coefficients / product(dx**derivatives(index_start:index_end))
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
      if(allocated(FDstencil_type_input%derivatives_order)) deallocate(FDstencil_type_input%derivatives_order)
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
      write(iounit, *) "derivatives_order: ", FDstencil_type_input%derivatives_order
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
      integer :: global_index, start_index, end_index, num_equations, num_taylor_elements, num_rhs, info
      integer, dimension(product(stencil_sizes)) :: ipiv
      real, dimension(product(stencil_sizes)) :: taylor_coefficients
      real, dimension((product(stencil_sizes)**(2*ndims))) :: coefficient_matrix

      info = 0

      num_equations = product(stencil_sizes)
      num_taylor_elements = num_equations
      num_rhs = 1

      ! Initialize the stencil coefficients. It starts as the rhs of the linear system.
      stencil_coefficients = 0.0
      stencil_coefficients(IDX_XD(ndims, stencil_sizes, derivative + 1)) = 1.0 ! derivative + 1 because the index should start at 1. We specify the derivative from 0.

      do global_index = 1, num_equations
         index = IDX_XD_INV(ndims, stencil_sizes, global_index)
         index = index - alphas - 1 ! go from alpha to beta
         start_index = (global_index - 1) * num_taylor_elements + 1
         end_index = global_index * num_taylor_elements

         ! The Taylor coefficients are calculated for each point of the stencil
         call calculate_taylor_coefficients(ndims, stencil_sizes, real(index,kind=8), taylor_coefficients)
         coefficient_matrix(start_index:end_index) = taylor_coefficients

      end do

      ! Call DGESV to solve the linear system
      call DGESV(num_equations, num_rhs, coefficient_matrix, num_taylor_elements, ipiv, &
         stencil_coefficients, num_taylor_elements, info)

      if (info /= 0) then
         print *, "DGESV reported an error: ", info
         stop
      end if

   end subroutine calculate_finite_difference_stencil

   !> Calculate the Taylor coefficients up to order s
   !! A_ij = sum_{i=0}^{p-1} sum_{j=0}^{q-1} 1/(i! * j!) (n*h)^i * (m*k)^j
   !! Here the factorial is slow. We could use the loop to calculate the factorial. But this is for another time. Should be fine for setup.
   pure subroutine calculate_taylor_coefficients(ndims, stencil_sizes, dx, taylor_coefficients)
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

         numerator = dx ** real(current_index-1,kind=8)
         denominator = gamma(real(current_index,kind=8))
         taylor_coefficients(global_index) = product(numerator / denominator)
      end do

   end subroutine calculate_taylor_coefficients

   !> Apply the finite difference stencil to a matrix
   pure function apply_FDstencil(ndims, FDstencil_type_input, dims, matrix, index) result(val)
      integer, intent(in) :: ndims
      type(FDstencil_type), intent(in) :: FDstencil_type_input
      integer, dimension(ndims), intent(in) :: dims
      real, dimension(product(dims)), intent(in) :: matrix
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

   !> subroutine to calculate the order of the derivative
   pure function calculate_derivative_order(ndims, num_derivatives, derivatives, stencil_sizes) result(order)
      integer, intent(in) :: ndims, num_derivatives
      integer, dimension(ndims*num_derivatives), intent(in) :: derivatives
      integer, dimension(ndims), intent(in) :: stencil_sizes

      integer, dimension(num_derivatives):: order

      integer :: global_index, index_start, index_end

      do global_index = 1, num_derivatives
         index_start = (global_index-1)*ndims+1
         index_end = global_index*ndims
         order(global_index) = minval(stencil_sizes - derivatives(index_start:index_end))
      end do

   end function calculate_derivative_order

   !> subroutine to find the requested order from the derivatives
   function find_stencil_size_from_order(ndims, num_derivatives, derivatives, order) result(stencil_sizes)
      integer, intent(in) :: ndims, num_derivatives, order
      integer, dimension(ndims*num_derivatives), intent(in) :: derivatives

      integer, dimension(ndims) :: stencil_sizes

      integer :: global_index, index_start, index_end

      stencil_sizes = 0

      do global_index = 1, num_derivatives
         index_start = (global_index-1)*ndims+1
         index_end = global_index*ndims
         stencil_sizes = max(stencil_sizes, derivatives(index_start:index_end) + order)
      end do

      stencil_sizes = stencil_sizes - 1

   end function find_stencil_size_from_order

end module FD_module
