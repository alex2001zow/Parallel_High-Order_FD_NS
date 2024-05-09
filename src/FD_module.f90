module FD_module
   use utility_functions_module, only : IDX_XD, IDX_XD_INV, print_real_array, sleeper_function
   use block_module, only : block_type
   implicit none

   private

   type FDstencil_type
      integer :: num_derivatives, num_stencil_elements, combination_of_stencil_sizes, stencil_coefficients_size
      integer, dimension(:), allocatable :: derivatives, stencil_sizes
      real, dimension(:), allocatable :: stencil_coefficients, scaled_stencil_coefficients
   end type FDstencil_type

   public :: FDstencil_type, create_finite_difference_stencils, deallocate_finite_difference_stencil, &
      print_finite_difference_stencil
   public :: apply_FDstencil, update_value_from_stencil
   public :: determine_alpha, alpha_2_global, global_2_start_end, get_FD_coefficients_from_index
   public :: calculate_scaled_coefficients

contains

   !> Create a finite difference stencil
   subroutine create_finite_difference_stencils(ndims, num_derivatives, derivatives, stencil_sizes, FDstencil_type_output)
      integer, intent(in) :: ndims, num_derivatives
      integer, dimension(ndims * num_derivatives), intent(in) :: derivatives
      integer, dimension(ndims), intent(in) :: stencil_sizes

      type(FDstencil_type), intent(out) :: FDstencil_type_output

      integer :: global_index, start_index, end_index
      integer, dimension(ndims) :: alphas, betas

      FDstencil_type_output%num_derivatives = num_derivatives
      FDstencil_type_output%num_stencil_elements = product(stencil_sizes)
      FDstencil_type_output%combination_of_stencil_sizes = product(stencil_sizes)
      FDstencil_type_output%stencil_coefficients_size = FDstencil_type_output%combination_of_stencil_sizes*&
         FDstencil_type_output%num_stencil_elements*num_derivatives

      allocate(FDstencil_type_output%derivatives(ndims * num_derivatives))
      allocate(FDstencil_type_output%stencil_sizes(ndims))
      allocate(FDstencil_type_output%stencil_coefficients(FDstencil_type_output%stencil_coefficients_size))
      allocate(FDstencil_type_output%scaled_stencil_coefficients(FDstencil_type_output%stencil_coefficients_size))

      FDstencil_type_output%derivatives = derivatives
      FDstencil_type_output%stencil_sizes = stencil_sizes

      do global_index = 1, FDstencil_type_output%combination_of_stencil_sizes
         call global_2_alpha_beta(ndims, FDstencil_type_output%stencil_sizes, global_index, alphas, betas)
         call global_2_start_end(ndims, num_derivatives, FDstencil_type_output%stencil_sizes, global_index, &
            start_index, end_index)

         call create_finite_difference_stencil_from_alpha_and_beta(ndims, num_derivatives, derivatives, alphas, betas, &
            FDstencil_type_output%stencil_coefficients(start_index:end_index))
      end do

   end subroutine create_finite_difference_stencils

   !> Create a finite difference stencil for a given value of alpha and beta
   subroutine create_finite_difference_stencil_from_alpha_and_beta(ndims, num_derivatives, derivatives, alphas, betas, &
      stencil_coefficients)
      integer, intent(in) :: ndims, num_derivatives
      integer, dimension(ndims * num_derivatives), intent(in) :: derivatives
      integer, dimension(ndims), intent(in) :: alphas, betas
      real, dimension(:), intent(out) :: stencil_coefficients

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

         stencil_coefficients_index_start = (global_index-1)*num_stencil_elements + 1
         stencil_coefficients_index_end = global_index*num_stencil_elements
         stencil_coefficients(stencil_coefficients_index_start:stencil_coefficients_index_end) = temp_stencil_coefficients
      end do

   end subroutine create_finite_difference_stencil_from_alpha_and_beta

   !> Destroy a finite difference stencil
   pure subroutine deallocate_finite_difference_stencil(FDstencil_type_input)
      type(FDstencil_type), intent(inout) :: FDstencil_type_input

      if(allocated(FDstencil_type_input%derivatives)) deallocate(FDstencil_type_input%derivatives)
      if(allocated(FDstencil_type_input%stencil_sizes)) deallocate(FDstencil_type_input%stencil_sizes)

      if(allocated(FDstencil_type_input%stencil_coefficients)) deallocate(FDstencil_type_input%stencil_coefficients)
      if(allocated(FDstencil_type_input%scaled_stencil_coefficients)) deallocate(FDstencil_type_input%scaled_stencil_coefficients)

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
   pure subroutine apply_FDstencil(ndims, num_elements, element_in_block, stencil_sizes, alphas, &
      coefficients, dims, index, matrix, val)
      integer, intent(in) :: ndims, num_elements, element_in_block
      integer, dimension(ndims), intent(in) :: stencil_sizes, alphas, dims, index
      real, dimension(product(stencil_sizes)), intent(in) :: coefficients
      real, dimension(product(dims) * num_elements), intent(in) :: matrix
      real, intent(out) :: val

      real :: stencil_val, matrix_val
      integer :: global_index, block_index

      integer, dimension(ndims) :: index_stencil

      val = 0.0

      ! Run through the stencil coefficients including the center coefficient
      do global_index = 1, product(stencil_sizes)
         call IDX_XD_INV(ndims, stencil_sizes, global_index, index_stencil)
         index_stencil = index_stencil - alphas - 1
         call IDX_XD(ndims, dims, index + index_stencil, block_index)
         block_index = block_index + element_in_block ! Shift the global index depending on the element in the block. Should be 0 for the first and 1 for the second and so on.

         stencil_val = coefficients(global_index)
         matrix_val = matrix(block_index)

         val = val + stencil_val * matrix_val
      end do

   end subroutine apply_FDstencil

   pure subroutine update_value_from_stencil(ndims, num_elements, element_in_block, stencil_sizes, alphas, &
      coefficients, dims, index, matrix, f_val, val)
      integer, intent(in) :: ndims, num_elements, element_in_block
      integer, dimension(ndims), intent(in) :: stencil_sizes, alphas, dims, index
      real, dimension(product(stencil_sizes)), intent(in) :: coefficients
      real, dimension(product(dims) * num_elements), intent(in) :: matrix
      real, intent(in) :: f_val
      real, intent(out) :: val

      integer :: global_matrix_index, center_coefficient_index
      real :: center_coefficient_value

      call IDX_XD(ndims, dims, index, global_matrix_index)
      call IDX_XD(ndims, stencil_sizes, alphas + 1, center_coefficient_index)

      center_coefficient_value = coefficients(center_coefficient_index)

      call apply_FDstencil(ndims, num_elements, element_in_block, stencil_sizes, alphas, coefficients, dims, index, matrix, val)

      val = val - center_coefficient_value * matrix(global_matrix_index + element_in_block) ! Subtract the center coefficient times the matrix value.
      val = (f_val - val) / center_coefficient_value ! Divide by the center coefficient to solve for the new value.

   end subroutine update_value_from_stencil

   !> Determine alpha and beta from a matrix index.
   pure subroutine determine_alpha(ndims, stencil_sizes, matrix_begin, matrix_end, matrix_index, alpha)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: stencil_sizes, matrix_begin, matrix_end, matrix_index
      integer, dimension(ndims), intent(out) :: alpha

      integer, dimension(ndims) :: elements_to_begin, elements_to_end, beta

      elements_to_begin = matrix_index - matrix_begin
      elements_to_end = matrix_end - matrix_index

      ! Try a centered stencil first
      alpha = min(elements_to_begin, stencil_sizes / 2)

      beta = (stencil_sizes - 1) - alpha
      alpha = alpha + max(beta - elements_to_end, 0)
      beta = (stencil_sizes - 1) - alpha ! Not needed but this will give the correct beta. We are really only interested in alpha.

   end subroutine determine_alpha

   pure subroutine global_2_alpha_beta(ndims, stencil_sizes, global_index, alphas, betas)
      integer, intent(in) :: ndims, global_index
      integer, dimension(ndims), intent(in) :: stencil_sizes
      integer, dimension(ndims), intent(out) :: alphas, betas

      integer, dimension(ndims) :: current_index

      call IDX_XD_INV(ndims, stencil_sizes, global_index, current_index)
      current_index = current_index - 1
      alphas = stencil_sizes - current_index - 1
      betas = current_index

   end subroutine global_2_alpha_beta

   pure subroutine alpha_2_global(ndims, stencil_sizes, alphas, global_index)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: stencil_sizes, alphas
      integer, intent(out) :: global_index

      call IDX_XD(ndims, stencil_sizes, stencil_sizes - alphas, global_index)

   end subroutine alpha_2_global

   pure subroutine global_2_start_end(ndims, num_derivatives, stencil_sizes, global_index, start_index, end_index)
      integer, intent(in) :: ndims, num_derivatives, global_index
      integer, dimension(ndims), intent(in) :: stencil_sizes
      integer, intent(out) :: start_index, end_index

      integer :: num_stencil_elements

      num_stencil_elements = product(stencil_sizes)

      start_index = (global_index-1)*num_stencil_elements*num_derivatives + 1
      end_index = global_index*num_stencil_elements*num_derivatives

   end subroutine global_2_start_end

   pure subroutine get_FD_coefficients_from_index(ndims, num_derivatives, stencil_sizes, start_dims, dims, index, &
      coefficients, alpha, pointer_to_coefficients)
      integer, intent(in) :: ndims, num_derivatives
      integer, dimension(ndims), intent(in) :: stencil_sizes, start_dims, dims, index
      real, dimension(:), target, intent(inout) :: coefficients
      integer, dimension(ndims), intent(out) :: alpha
      real, dimension(:), pointer, intent(out) :: pointer_to_coefficients

      integer :: global_index, start_index, end_index

      ! Determine the alpha for the current index
      call determine_alpha(ndims, stencil_sizes, start_dims, dims, index, alpha)

      ! Find the coefficients from alpha.
      call alpha_2_global(ndims, stencil_sizes, alpha, global_index)
      call global_2_start_end(ndims, num_derivatives, stencil_sizes, global_index, start_index, end_index)

      pointer_to_coefficients => coefficients(start_index:end_index)

   end subroutine get_FD_coefficients_from_index

   !> Routine to scale the finite difference coefficients depending on the grid spacing(dx)
   pure subroutine calculate_scaled_coefficients(ndims, dx, FDstencil_type_input)
      integer, intent(in) :: ndims
      real, dimension(ndims), intent(in) :: dx
      type(FDstencil_type), target, intent(inout) :: FDstencil_type_input

      integer :: global_index, start_index, end_index
      integer :: derivative_global_index, derivative_start_index, derivative_end_index
      integer :: coefficient_start_index, coefficient_end_index
      integer, dimension(ndims) :: alphas, betas, derivative
      integer :: num_stencil_elements

      real :: scale

      real, dimension(:), pointer :: ptr_org, ptr_scaled

      ptr_org => FDstencil_type_input%stencil_coefficients
      ptr_scaled => FDstencil_type_input%scaled_stencil_coefficients

      num_stencil_elements = product(FDstencil_type_input%stencil_sizes)

      do global_index = 1, FDstencil_type_input%combination_of_stencil_sizes
         call global_2_alpha_beta(ndims, FDstencil_type_input%stencil_sizes, global_index, alphas, betas)
         call global_2_start_end(ndims, FDstencil_type_input%num_derivatives, FDstencil_type_input%stencil_sizes, global_index, &
            start_index, end_index)

         do derivative_global_index = 1, FDstencil_type_input%num_derivatives
            derivative_start_index = (derivative_global_index-1)*ndims + 1
            derivative_end_index = derivative_global_index*ndims
            derivative = FDstencil_type_input%derivatives(derivative_start_index:derivative_end_index)

            scale = product(dx**derivative)

            coefficient_start_index = (derivative_global_index-1)*num_stencil_elements
            coefficient_end_index = derivative_global_index*num_stencil_elements - 1

            ptr_scaled(start_index+coefficient_start_index:start_index+coefficient_end_index) = &
               ptr_org(start_index+coefficient_start_index:start_index+coefficient_end_index) / scale
         end do
      end do

   end subroutine calculate_scaled_coefficients

end module FD_module
