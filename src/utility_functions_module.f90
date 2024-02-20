module utility_functions_module
   implicit none

   private
   public :: IDX_XD, IDX_XD_INV, get_indices, sleeper_function

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

   !> Local indices from ndims, dims, global_index.
   !! 1D: ndims = 1, dims = [N], global_index = i, indices = [i]
   !! 2D: ndims = 2, dims = [N, M], global_index = i + N * (j - 1), indices = [i, j]
   !! 3D: ndims = 3, dims = [N, M, K], global_index = i + N * (j - 1) + N * M * (k - 1), indices = [i, j, k]
   !! 4D: ndims = 4, dims = [N, M, K, L], global_index = i + N * (j - 1) + N * M * (k - 1) + N * M * K * (l - 1), indices = [i, j, k, l]
   !! etc.
   function IDX_XD_INV(ndims, dims, global_index_in) result(indices)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims
      integer, intent(in) :: global_index_in
      integer, dimension(ndims) :: indices
      integer :: d, product, global_index

      ! Adjust index for 1-based indexing of Fortran
      global_index = global_index_in - 1

      product = 1
      do d = 1, ndims - 1
         product = product * dims(d)
      end do

      do d = ndims, 2, -1
         indices(d) = global_index / product + 1
         global_index = mod(global_index, product)
         product = product / dims(d - 1)
      end do

      indices(1) = global_index + 1
   end function IDX_XD_INV

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
