module utility_functions_module
   implicit none

   private
   public :: IDX_XD, IDX_XD_REVERSED, get_indices, sleeper_function, print_matrix, flip_1d_integer_array

contains

   !> Global index from ndims, dims, indices.
   !! 1D: ndims = 1, dims = [N], indices = [i], global_index = i
   !! 2D: ndims = 2, dims = [N, M], indices = [i, j], global_index = i + N * (j - 1)
   !! 3D: ndims = 3, dims = [N, M, K], indices = [i, j, k], global_index = i + N * (j - 1) + N * M * (k - 1)
   !! 4D: ndims = 4, dims = [N, M, K, L], indices = [i, j, k, l], global_index = i + N * (j - 1) + N * M * (k - 1) + N * M * K * (l - 1)
   !! etc.
   function IDX_XD(ndims, dims, indices) result(global_index)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims, indices
      integer :: global_index
      integer :: d, product

      global_index = indices(1) - 1
      product = 1

      do d = 2, ndims
         product = product * dims(d - 1)
         global_index = global_index + (indices(d) - 1) * product
      end do

      ! Adjust index for 1-based indexing of Fortran
      global_index = global_index + 1
   end function IDX_XD

   !> Global index from ndims, dims, indices. Reversed order. More intuitive for 2D and 3D arrays. I think.
   !! Everything using IDX_XD should be replaced with IDX_XD_REVERSED eventually.
   function IDX_XD_REVERSED(ndims, dims, indices) result(global_index)
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
   end function IDX_XD_REVERSED


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
               indices(global_index) = IDX_XD(ndims, dims, [jj, ii])
               global_index = global_index + 1
            end do
         end do
      end if

      if(ndims == 3) then
         global_index = 1
         do ii = begin(1), end(1)
            do jj = begin(2),end(2)
               do kk = begin(3),end(2)
                  indices(global_index) = IDX_XD(ndims, dims, [kk, jj, ii])
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
               global_index = IDX_XD(ndims, dims, [jj, ii])
               write(*, '(F10.3, " ")', advance="no") matrix(global_index)
            end do
            write(*, *)
         end do
      end if

   end subroutine print_matrix

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

   !> Routine to flip a 1D array.
   subroutine flip_1d_integer_array(array_size, array)
      integer, intent(in) :: array_size
      integer, dimension(array_size), intent(inout) :: array

      integer :: temp
      integer :: ii

      do ii = 1, array_size / 2
         temp = array(ii)
         array(ii) = array(array_size - ii + 1)
         array(array_size - ii + 1) = temp
      end do
   end subroutine flip_1d_integer_array

end module utility_functions_module
