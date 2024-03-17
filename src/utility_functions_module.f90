module utility_functions_module
   implicit none

   private
   public :: IDX_XD, IDX_XD_INV, print_real_array, print_integer_array, read_input_from_command_line, &
      swap_pointers, sleeper_function

contains

   !> Global index from ndims, dims, indices.
   pure function IDX_XD(ndims, dims, indices) result(global_index)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: dims, indices
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
      integer, dimension(:), intent(in) :: dims
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

   !> Routine to print an array of type real
   subroutine print_real_array(ndims, dims, array, stride, title, iounit)
      integer, intent(in) :: ndims, stride, iounit
      integer, dimension(:), intent(in) :: dims
      real, dimension(:), intent(in) :: array
      character(len=*), intent(in) :: title

      integer :: global_index, d, start_index, end_index
      integer, dimension(ndims) :: indices
      real, dimension(stride) :: temp_array
      character(len=255) :: format_string

      write(format_string, '("(", i0, "F", i0, ".", i0, ")")') stride, kind(array(1)), 3

      write(iounit, *) title
      do global_index = 1, product(dims)
         indices = IDX_XD_INV(ndims, dims, global_index)
         start_index = (global_index - 1) * stride + 1
         end_index = global_index * stride
         temp_array = array(start_index:end_index)
         write(iounit, format_string, advance='no') temp_array

         do d = ndims, 2, -1
            if (indices(d) == dims(d)) then
               write(iounit, *)
               exit
            end if
         end do
      end do

   end subroutine print_real_array

   !> Routine to print an array of type integer
   subroutine print_integer_array(ndims, dims, array, stride, title, iounit)
      integer, intent(in) :: ndims, stride, iounit
      integer, dimension(:), intent(in) :: dims
      integer, dimension(:), intent(in) :: array
      character(len=*), intent(in) :: title

      integer :: global_index, d, start_index, end_index
      integer, dimension(ndims) :: indices
      integer, dimension(stride) :: temp_array
      character(len=255) :: format_string

      write(format_string, '("(", i0, "I", i0, ")")') stride, kind(array(1))

      write(iounit, *) title
      do global_index = 1, product(dims)
         indices = IDX_XD_INV(ndims, dims, global_index)
         start_index = (global_index - 1) * stride + 1
         end_index = global_index * stride
         temp_array = array(start_index:end_index)
         write(iounit, format_string, advance='no') temp_array

         do d = ndims, 2, -1
            if (indices(d) == dims(d)) then
               write(iounit, *)
               exit
            end if
         end do
      end do

   end subroutine print_integer_array

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
