module utility_functions_module
   implicit none

   private
   public :: IDX_XD, IDX_XD_INV
   public :: print_real_array, print_integer_array
   public :: open_txt_file, close_txt_file, read_input_from_command_line
   public :: calculate_dx, swap_pointers
   public :: sleeper_function

contains

   !> Global index from ndims, dims, indices.
   pure subroutine IDX_XD(ndims, dims, indices, global_index)
      integer, intent(in) :: ndims
      integer, dimension(:), intent(in) :: dims, indices
      integer, intent(out) :: global_index

      integer :: d, stride

      global_index = indices(ndims) - 1
      stride = 1

      do d = ndims-1,1,-1
         stride = stride * dims(d + 1)
         global_index = global_index + (indices(d) - 1) * stride
      end do

      ! Adjust index for 1-based indexing of Fortran
      global_index = global_index + 1
   end subroutine IDX_XD

   !> Dimensional indices from ndims, dims, global_index.
   pure subroutine IDX_XD_INV(ndims, dims, global_index, indices)
      integer, intent(in) :: ndims, global_index
      integer, dimension(:), intent(in) :: dims
      integer, dimension(:), intent(out) :: indices

      integer :: d, remaining_index, stride

      ! Adjust index for 1-based indexing of Fortran
      remaining_index = global_index - 1

      ! Start with the product of all dimensions except the last one
      stride = product(dims(2:ndims))
      !stride = 1
      !do d = 2, ndims
      !   stride = stride * dims(d)
      !end do

      do d = 1, ndims-1
         indices(d) = (remaining_index / stride) + 1
         remaining_index = mod(remaining_index, stride)
         stride = stride / dims(d+1)
      end do

      ! Handle the last dimension separately
      indices(ndims) = remaining_index + 1

   end subroutine IDX_XD_INV

   !> Routine to print an array of type real
   subroutine print_real_array(ndims, dims, array, stride, title, iounit)
      integer, intent(in) :: ndims, stride, iounit
      integer, dimension(:), intent(in) :: dims
      real, dimension(:), intent(in) :: array
      character(len=*), intent(in) :: title

      integer :: do_end, global_index, d, start_index, end_index
      integer, dimension(ndims) :: indices

      character(len=32) :: format_str   ! Length to be adjusted based on expected size
      character(len=10) :: element_format
      character(len=255) :: write_format

      element_format = 'F10.3'  ! Format for each element
      write(format_str, '(I32)') dims(1) * stride
      write(write_format, '(*(A))') '(', trim(adjustl(format_str)), '(', element_format, '))'


      if(ndims > 1) then
         do_end = product(dims(2:ndims))
      else
         do_end = 1
      end if

      write(iounit, *) title
      do global_index = 1, do_end
         call IDX_XD_INV(ndims, dims, global_index, indices)
         start_index = (global_index - 1) * dims(1) * stride + 1
         end_index = global_index * dims(1) * stride
         write(iounit, write_format) array(start_index:end_index)

         do d = ndims, 3, -1
            if (indices(d) == dims(d)) then
               write(iounit, *)
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

      integer :: do_end, global_index, d, start_index, end_index
      integer, dimension(ndims) :: indices

      character(len=32) :: format_str    ! Length to be adjusted based on expected size
      character(len=5) :: element_format
      character(len=255) :: write_format

      element_format = 'I6'  ! Format for each element, adjust width as needed
      write(format_str, '(I32)') dims(1) * stride
      write(write_format, '(*(A))') '(', trim(adjustl(format_str)), '(', element_format, '))'


      if(ndims > 1) then
         do_end = product(dims(2:ndims))
      else
         do_end = 1
      end if

      write(iounit, *) title
      do global_index = 1, do_end
         call IDX_XD_INV(ndims, dims, global_index, indices)
         start_index = (global_index - 1) * dims(1) * stride + 1
         end_index = global_index * dims(1) * stride
         write(iounit, write_format) array(start_index:end_index)

         do d = ndims, 3, -1
            if (indices(d) == dims(d)) then
               write(iounit, *)
            end if
         end do
      end do

   end subroutine print_integer_array

   !> Routine to open a file for writing
   subroutine open_txt_file(filename, rank, iounit)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: rank

      integer :: iounit, ios
      character(255) :: file_with_rank

      ! Create a modified filename by appending the rank to the original filename
      write(file_with_rank, '(A, I0, A)') trim(filename), rank, ".txt"

      ! Open the file for writing, associate it with a logical unit (iounit)
      open(newunit=iounit, file=file_with_rank, status='replace', action='write', iostat=ios)

      ! Check for errors in opening the file
      if (ios /= 0) then
         print *, "Error opening file: ", file_with_rank
         return
      end if

   end subroutine open_txt_file

   !> Subroutine to close a file
   subroutine close_txt_file(iounit)
      integer, intent(in) :: iounit

      close(iounit)
   end subroutine close_txt_file

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

   !> Routine to calculate the grid spacing dx. This can return negative values if the domain is reversed. Pretty sure we always want positive values due to the way we calculate the stencils.
   pure subroutine calculate_dx(domain_begin, domain_end, grid_size_with_ghosts, dx)
      real, dimension(:), intent(in) :: domain_begin, domain_end
      integer, dimension(:), intent(in) :: grid_size_with_ghosts
      real, dimension(:), intent(out) :: dx

      dx = (domain_end - domain_begin) / (grid_size_with_ghosts - 1)
   end subroutine calculate_dx

   !> Subroutine to swap two pointers
   pure subroutine swap_pointers(ptr1, ptr2)
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
