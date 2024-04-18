module FD_module
   use utility_functions_module, only : IDX_XD, IDX_XD_INV, print_real_array, sleeper_function
   use block_module, only : block_type
   implicit none

   private

   type FDstencil_type
      integer :: num_derivatives, num_stencil_elements, combination_of_stencil_sizes, stencil_coefficients_size
      integer, dimension(:), allocatable :: derivatives, derivatives_order, stencil_sizes
      real, dimension(:), allocatable :: dx, stencil_coefficients
   end type FDstencil_type

   public :: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil
   public :: apply_FDstencil, determine_alpha, alpha_2_global, global_2_start_end

contains

   !> Create a finite difference stencil
   subroutine create_finite_difference_stencils(ndims, num_derivatives, derivatives, dx, stencil_sizes, FDstencil_type_output)
      integer, intent(in) :: ndims, num_derivatives
      integer, dimension(ndims * num_derivatives), intent(in) :: derivatives
      integer, dimension(ndims), intent(in) :: stencil_sizes
      real,  dimension(ndims), intent(in) :: dx

      type(FDstencil_type), intent(out) :: FDstencil_type_output

      integer :: global_index, start_index, end_index
      integer, dimension(ndims) :: alphas, betas

      FDstencil_type_output%num_derivatives = num_derivatives
      FDstencil_type_output%num_stencil_elements = product(stencil_sizes)
      FDstencil_type_output%combination_of_stencil_sizes = product(stencil_sizes)
      FDstencil_type_output%stencil_coefficients_size = FDstencil_type_output%combination_of_stencil_sizes*&
         FDstencil_type_output%num_stencil_elements*num_derivatives

      allocate(FDstencil_type_output%derivatives(ndims * num_derivatives))
      allocate(FDstencil_type_output%derivatives_order(num_derivatives))
      allocate(FDstencil_type_output%dx(ndims))
      allocate(FDstencil_type_output%stencil_sizes(ndims))
      allocate(FDstencil_type_output%stencil_coefficients(FDstencil_type_output%stencil_coefficients_size))

      FDstencil_type_output%derivatives = derivatives
      FDstencil_type_output%dx = dx
      FDstencil_type_output%stencil_sizes = stencil_sizes
      call calculate_derivative_order(ndims, num_derivatives, derivatives, &
         FDstencil_type_output%stencil_sizes, FDstencil_type_output%derivatives_order)

      do global_index = 1, FDstencil_type_output%combination_of_stencil_sizes
         call global_2_alpha_beta(ndims, FDstencil_type_output%stencil_sizes, global_index, alphas, betas)

         call global_2_start_end(ndims, num_derivatives, FDstencil_type_output%stencil_sizes, global_index, &
            start_index, end_index)

         call create_finite_difference_stencil_from_alpha_and_beta(ndims, num_derivatives, derivatives, alphas, betas, dx, &
            FDstencil_type_output%stencil_coefficients(start_index:end_index))
      end do

   end subroutine create_finite_difference_stencils

   !> Create a finite difference stencil for a given value of alpha and beta
   subroutine create_finite_difference_stencil_from_alpha_and_beta(ndims, num_derivatives, derivatives, alphas, betas, dx, &
      stencil_coefficients)
      integer, intent(in) :: ndims, num_derivatives
      integer, dimension(ndims * num_derivatives), intent(in) :: derivatives
      integer, dimension(ndims), intent(in) :: alphas, betas
      real, dimension(ndims), intent(in) :: dx
      real, dimension(:), intent(inout) :: stencil_coefficients

      integer :: global_index, index_start, index_end
      integer :: stencil_coefficients_index_start, stencil_coefficients_index_end
      integer, dimension(ndims) :: stencil_sizes, derivative
      integer :: num_stencil_elements
      real, dimension(product(alphas + 1 + betas)) :: temp_stencil_coefficients

      stencil_sizes = alphas + 1 + betas
      num_stencil_elements = product(stencil_sizes)

      do global_index = 1, num_derivatives
         index_start = (global_index-1)*ndims+1
         index_end = global_index*ndims
         derivative = derivatives(index_start:index_end)
         call calculate_finite_difference_stencil(ndims, alphas, betas, derivative, temp_stencil_coefficients)
         !temp_stencil_coefficients = temp_stencil_coefficients / product(dx**derivative) ! This should be removed and calculates on a per block basis. Then we can use multigrid.

         stencil_coefficients_index_start = (global_index-1)*num_stencil_elements + 1
         stencil_coefficients_index_end = global_index*num_stencil_elements
         stencil_coefficients(stencil_coefficients_index_start:stencil_coefficients_index_end) &
            = temp_stencil_coefficients
      end do

   end subroutine create_finite_difference_stencil_from_alpha_and_beta

   !> Destroy a finite difference stencil
   pure subroutine deallocate_finite_difference_stencil(FDstencil_type_input)
      type(FDstencil_type), intent(inout) :: FDstencil_type_input

      if(allocated(FDstencil_type_input%derivatives)) deallocate(FDstencil_type_input%derivatives)
      if(allocated(FDstencil_type_input%derivatives_order)) deallocate(FDstencil_type_input%derivatives_order)
      if(allocated(FDstencil_type_input%dx)) deallocate(FDstencil_type_input%dx)
      if(allocated(FDstencil_type_input%stencil_sizes)) deallocate(FDstencil_type_input%stencil_sizes)

      if(allocated(FDstencil_type_input%stencil_coefficients)) deallocate(FDstencil_type_input%stencil_coefficients)

   end subroutine deallocate_finite_difference_stencil

   !> Print the finite difference stencil
   subroutine print_finite_difference_stencil(ndims, FDstencil_type_input, iounit)
      integer, intent(in) :: ndims, iounit
      type(FDstencil_type), intent(in) :: FDstencil_type_input

      integer :: combination_of_stencil_sizes, num_stencil_elements, num_derivatives
      integer :: alpha_beta_global_index, alpha_beta_start_index, alpha_beta_end_index
      integer :: derivative_global_index, derivative_start_index, derivative_end_index
      integer, dimension(ndims) :: current_index, alphas, betas

      combination_of_stencil_sizes = FDstencil_type_input%combination_of_stencil_sizes
      num_stencil_elements = FDstencil_type_input%num_stencil_elements
      num_derivatives = FDstencil_type_input%num_derivatives

      ! Newlines for readability
      write(iounit, *)
      write(iounit, *) "FDstencil_type: "
      write(iounit, *)

      write(iounit, *) "num_derivatives: ", FDstencil_type_input%num_derivatives
      write(iounit, *) "derivatives: ", FDstencil_type_input%derivatives
      write(iounit, *) "derivatives_order: ", FDstencil_type_input%derivatives_order
      write(iounit, *) "dx: ", FDstencil_type_input%dx
      write(iounit, *) "stencil_sizes: ", FDstencil_type_input%stencil_sizes

      write(iounit, *) "stencil_coefficients: "
      do alpha_beta_global_index = 1, combination_of_stencil_sizes
         call IDX_XD_INV(ndims, FDstencil_type_input%stencil_sizes, alpha_beta_global_index, current_index)
         current_index = current_index - 1
         alphas = FDstencil_type_input%stencil_sizes - current_index - 1
         betas = current_index

         alpha_beta_start_index = (alpha_beta_global_index-1)*num_stencil_elements*num_derivatives + 1
         alpha_beta_end_index = alpha_beta_global_index*num_stencil_elements*num_derivatives

         write(iounit, *) "Alpha and Beta: ", alphas, betas
         do derivative_global_index = 1, num_derivatives

            derivative_start_index = (derivative_global_index-1)*num_stencil_elements
            derivative_end_index = derivative_global_index*num_stencil_elements-1

            write(iounit, *) "Derivative: ", FDstencil_type_input%derivatives(&
               (derivative_global_index-1)*ndims+1:derivative_global_index*ndims)
            call print_real_array(ndims, FDstencil_type_input%stencil_sizes, &
               FDstencil_type_input%stencil_coefficients(&
               (alpha_beta_start_index+derivative_start_index):(alpha_beta_start_index+derivative_end_index)), 1, "", iounit)
         end do
      end do

   end subroutine print_finite_difference_stencil

   !> Calculate the finite difference stencil coefficients.
   !! sum_{n=-alpha_1}^{beta_1} sum_{m=-alpha_2}^{beta_2} a_nm * A_ij
   subroutine calculate_finite_difference_stencil(ndims, alphas, betas, derivative, stencil_coefficients)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: alphas, betas, derivative

      integer, dimension(ndims) :: stencil_sizes

      real, dimension(product(alphas + 1 + betas)), intent(out) :: stencil_coefficients

      integer, dimension(ndims) :: index
      integer :: global_index, start_index, end_index, num_equations, num_taylor_elements, num_rhs, info
      integer, dimension(product(alphas + 1 + betas)) :: ipiv
      real, dimension(product(alphas + 1 + betas)) :: taylor_coefficients
      real, dimension(product(alphas + 1 + betas)**2) :: coefficient_matrix

      stencil_sizes = alphas + 1 + betas

      info = 0

      num_equations = product(stencil_sizes)
      num_taylor_elements = num_equations
      num_rhs = 1

      ! Initialize the stencil coefficients. It starts as the rhs of the linear system.
      call IDX_XD(ndims, stencil_sizes, derivative + 1, global_index) ! derivative + 1 because the index should start at 1. We specify the derivative from 0.
      stencil_coefficients = 0.0
      stencil_coefficients(global_index) = 1.0

      do global_index = 1, num_equations
         call IDX_XD_INV(ndims, stencil_sizes, global_index, index)
         index = index - alphas - 1 ! go from alpha to beta
         start_index = (global_index - 1) * num_taylor_elements + 1
         end_index = global_index * num_taylor_elements

         ! The Taylor coefficients are calculated for each point of the stencil
         call calculate_taylor_coefficients(ndims, stencil_sizes, real(index,kind=kind(taylor_coefficients(1))), &
            taylor_coefficients)
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
         call IDX_XD_INV(ndims, stencil_sizes, global_index, current_index)

         numerator = dx ** real(current_index-1,kind=8)
         denominator = gamma(real(current_index,kind=8))
         taylor_coefficients(global_index) = product(numerator / denominator)
      end do

   end subroutine calculate_taylor_coefficients

   !> Apply the finite difference stencil to a matrix
   pure subroutine apply_FDstencil(ndims, stencil_sizes, alphas, coefficients, dims, index, matrix, val)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: stencil_sizes, alphas, dims, index
      real, dimension(product(stencil_sizes)), intent(in) :: coefficients
      real, dimension(product(dims)), intent(in) :: matrix
      real, intent(inout) :: val

      real :: stencil_val, matrix_val
      integer :: global_index, block_index

      integer, dimension(ndims) :: index_stencil

      val = 0.0

      ! Run through the stencil coefficients including the center coefficient
      do global_index = 1, product(stencil_sizes)
         call IDX_XD_INV(ndims, stencil_sizes, global_index, index_stencil)
         index_stencil = index_stencil - alphas - 1
         call IDX_XD(ndims, dims, index + index_stencil, block_index)

         stencil_val = coefficients(global_index)
         matrix_val = matrix(block_index)

         val = val + stencil_val * matrix_val
      end do

   end subroutine apply_FDstencil

   !> DOUBLE CHECK THIS FUNCTION! Determine alpha and beta from a matrix index.
   pure subroutine determine_alpha(ndims, stencil_sizes, matrix_begin, matrix_end, matrix_index, alpha)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: stencil_sizes, matrix_begin, matrix_end, matrix_index
      integer, dimension(ndims), intent(inout) :: alpha

      integer, dimension(ndims) :: elements_to_begin, elements_to_end, centered_stencil, beta

      elements_to_begin = matrix_index - matrix_begin
      elements_to_end = matrix_end - matrix_index

      ! Try a centered stencil first
      centered_stencil = stencil_sizes / 2

      alpha = centered_stencil
      alpha = min(elements_to_begin, alpha)

      beta = (stencil_sizes - 1) - alpha
      alpha = alpha + max(beta - elements_to_end, 0)

   end subroutine determine_alpha

   !> Subroutine to calculate the order of the derivative
   pure subroutine calculate_derivative_order(ndims, num_derivatives, derivatives, stencil_sizes, order)
      integer, intent(in) :: ndims, num_derivatives
      integer, dimension(ndims*num_derivatives), intent(in) :: derivatives
      integer, dimension(ndims), intent(in) :: stencil_sizes

      integer, dimension(num_derivatives), intent(inout) :: order

      integer :: global_index, index_start, index_end

      do global_index = 1, num_derivatives
         index_start = (global_index-1)*ndims+1
         index_end = global_index*ndims
         order(global_index) = minval(stencil_sizes - derivatives(index_start:index_end))
      end do

   end subroutine calculate_derivative_order

   pure subroutine global_2_alpha_beta(ndims, stencil_sizes, global_index, alphas, betas)
      integer, intent(in) :: ndims, global_index
      integer, dimension(ndims), intent(in) :: stencil_sizes

      integer, dimension(ndims), intent(inout) :: alphas, betas

      integer, dimension(ndims) :: current_index

      call IDX_XD_INV(ndims, stencil_sizes, global_index, current_index)
      current_index = current_index - 1
      alphas = stencil_sizes - current_index - 1
      betas = current_index

   end subroutine global_2_alpha_beta

   pure subroutine alpha_2_global(ndims, stencil_sizes, alphas, global_index)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: stencil_sizes, alphas

      integer, intent(inout) :: global_index

      call IDX_XD(ndims, stencil_sizes, stencil_sizes - alphas, global_index)

   end subroutine alpha_2_global

   pure subroutine global_2_start_end(ndims, num_derivatives, stencil_sizes, global_index, start_index, end_index)
      integer, intent(in) :: ndims, num_derivatives, global_index
      integer, dimension(ndims), intent(in) :: stencil_sizes

      integer, intent(inout) :: start_index, end_index

      integer :: num_stencil_elements

      num_stencil_elements = product(stencil_sizes)

      start_index = (global_index-1)*num_stencil_elements*num_derivatives + 1
      end_index = global_index*num_stencil_elements*num_derivatives

   end subroutine global_2_start_end

end module FD_module
