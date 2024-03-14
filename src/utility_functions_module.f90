module utility_functions_module
   implicit none

   private
   public :: IDX_XD, IDX_XD_INV, print_matrix, read_input_from_command_line, find_abs_diff_matrices, swap_pointers, sleeper_function

contains

   !> Global index from ndims, dims, indices.
   pure function IDX_XD(ndims, dims, indices) result(global_index)
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

   !> Dimensional indices from ndims, dims, global_index.
   pure function IDX_XD_INV(ndims, dims, global_index) result(indices)
      integer, intent(in) :: ndims
      integer, dimension(ndims), intent(in) :: dims
      integer, intent(in) :: global_index

      integer, dimension(ndims) :: indices
      integer :: d, remaining_index, product

      ! Adjust index for 1-based indexing of Fortran
      remaining_index = global_index - 1

      ! Start with the product of all dimensions except the last one
      product = 1
      do d = 2, ndims
         product = product * dims(d)
      end do

      do d = 1, ndims-1
         indices(d) = (remaining_index / product) + 1
         remaining_index = mod(remaining_index, product)
         product = product / dims(d+1)
      end do

      ! Handle the last dimension separately
      indices(ndims) = remaining_index + 1

   end function IDX_XD_INV

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

   !> Routine to read input from the command line
   subroutine read_input_from_command_line(ndims, grid_size, processor_dim)
      integer, intent(out) :: ndims
      integer, dimension(:), allocatable, intent(out) :: grid_size, processor_dim

      character(255) :: temp_arg
      integer :: ii

      call get_command_argument(1, temp_arg)
      read(temp_arg,"(I5)") ndims

      allocate(grid_size(ndims))
      allocate(processor_dim(ndims))

      do ii = 1, ndims
         call get_command_argument(ii + 1, temp_arg)
         read(temp_arg,"(I5)") grid_size(ii)

         call get_command_argument(ii + (ndims + 1), temp_arg)
         read(temp_arg,"(I5)") processor_dim(ii)

         ! We have make sure that the inputs are properly divisible
         if(mod(grid_size(ii), processor_dim(ii)) /= 0) then
            print *, "Grid size is not divisible by the number of processors in dimension ", ii
            stop
         end if
      end do
   end subroutine read_input_from_command_line

   subroutine swap_pointers(ptr1, ptr2)
      real, dimension(:), pointer, intent(inout) :: ptr1, ptr2
      real, dimension(:), pointer :: temp_ptr

      temp_ptr => ptr1
      ptr1 => ptr2
      ptr2 => temp_ptr
   end subroutine swap_pointers

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
