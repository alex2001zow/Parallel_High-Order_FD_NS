module utility_functions_module
   implicit none

   private
   public :: IDX_XD, get_indices, print_matrix, sleeper_function, find_abs_diff_matrices

contains

   !> Global index from ndims, dims, indices.
   function IDX_XD(ndims, dims, indices) result(global_index)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims, indices
      integer :: global_index
      integer :: d, product

      global_index = indices(ndims) - 1
      product = 1

      do d = ndims-1, 1, -1
         product = product * dims(d + 1)
         global_index = global_index + (indices(d) - 1) * product
      end do

      ! Adjust index for 1-based indexing of Fortran
      global_index = global_index + 1
   end function IDX_XD


   !> Routine to get the indices for certain loop values OBS! MAKE SURE IT IS CORRECT. I THINK WE DO IT CORRECTLY BUT I AM NOT SURE
   subroutine get_indices(ndims, dims, begin, end, indices)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims, begin, end
      integer, allocatable, intent(out) :: indices(:)
      integer :: ii, jj, kk, global_index, number_of_indices

      number_of_indices = 1
      do ii = 1,ndims
         number_of_indices = number_of_indices * (end(ii) - begin(ii) + 1)
      end do

      allocate(indices(number_of_indices))

      global_index = 1

      if(ndims == 2) then
         do ii = begin(1), end(1)
            do jj = begin(2),end(2)
               indices(global_index) = IDX_XD(ndims, dims, [ii, jj])
               global_index = global_index + 1
            end do
         end do
      end if

      if(ndims == 3) then
         global_index = 1
         do ii = begin(1), end(1)
            do jj = begin(2),end(2)
               do kk = begin(3),end(2)
                  indices(global_index) = IDX_XD(ndims, dims, [ii, jj, kk])
                  global_index = global_index + 1
               end do
            end do
         end do
      end if

   end subroutine get_indices

   !> Routine to print a matrix.
   subroutine print_matrix(ndims, dims, matrix)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims
      real, dimension(product(dims)), intent(in) :: matrix

      integer :: ii, jj, global_index

      write(*, *) "size of matrix: ", size(matrix), "dims: ", dims, "product of dims", product(dims)

      if(ndims == 2) then
         do ii = 1, dims(1)
            do jj = 1, dims(2)
               global_index = IDX_XD(ndims, dims, [ii, jj])
               write(*, '(F10.3, " ")', advance="no") matrix(global_index)
            end do
            write(*, *)
         end do
      end if

   end subroutine print_matrix

   !> Routine to print the absolute difference between two matrices.
   subroutine find_abs_diff_matrices(ndims, dims, matrix1, matrix2)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims
      real, dimension(product(dims)), intent(in) :: matrix1, matrix2

      integer :: ii, jj, global_index

      if(ndims == 2) then
         do ii = 1, dims(1)
            do jj = 1, dims(2)
               global_index = IDX_XD(ndims, dims, [ii, jj])
               write(*, '(F10.3, " ")', advance="no") abs(matrix1(global_index) - matrix2(global_index))
            end do
            write(*, *)
         end do
      end if

   end subroutine find_abs_diff_matrices

   !> Routine to sleep for a certain amount of time. Used for debugging purposes.
   subroutine sleeper_function(sleep_time)
      implicit none
      integer, intent(in) :: sleep_time
      integer :: sleeper

      sleeper = 0
      do while(sleeper == 0)
         call sleep(sleep_time)
      end do
   end subroutine sleeper_function

end module utility_functions_module
