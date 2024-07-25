module utility_functions_module
   use, intrinsic :: iso_c_binding
   implicit none

   private
   public :: reshape_real_1D_to_2D, reshape_real_1D_to_3D, reshape_real_1D_to_4D, reshape_real_1D_to_5D, &
      reshape_real_1D_to_6D, reshape_real_1D_to_7D
   public :: reshape_integer_1D_to_2D, reshape_integer_1D_to_3D
   public :: IDX_XD, IDX_XD_INV
   public :: print_real_array, print_real_1D_array, print_real_2D_array, print_real_3D_array
   public :: print_integer_array, print_integer_1D_array, print_integer_2D_array
   public :: open_txt_file, close_txt_file, read_input_from_command_line
   public :: calculate_dx, calculate_CFL, calculate_dt_from_CFL, swap_pointers, swap_pointers_2D
   public :: sleeper_function

contains

   !> Routine to use c_f_pointer to reshape a contiguous real 1D-array into an 2D-array
   subroutine reshape_real_1D_to_2D(dims, array_1D, array_2D)
      integer, dimension(:), intent(in) :: dims
      real, contiguous, dimension(:), pointer, intent(in) :: array_1D
      real, contiguous, dimension(:,:), pointer, intent(out) :: array_2D

      ! Assign the C pointer to the target Fortran array
      call c_f_pointer(c_loc(array_1D), array_2D, dims)
   end subroutine reshape_real_1D_to_2D

   !> Routine to use c_f_pointer to reshape a contiguous real 1D-array into an 3D-array
   subroutine reshape_real_1D_to_3D(dims, array_1D, array_3D)
      integer, dimension(:), intent(in) :: dims
      real, contiguous, dimension(:), pointer, intent(in) :: array_1D
      real, contiguous, dimension(:,:,:), pointer, intent(out) :: array_3D

      ! Assign the C pointer to the target Fortran array
      call c_f_pointer(c_loc(array_1D), array_3D, dims)
   end subroutine reshape_real_1D_to_3D

   !> Routine to use c_f_pointer to reshape a contiguous real 1D-array into an 4D-array
   subroutine reshape_real_1D_to_4D(dims, array_1D, array_4D)
      integer, dimension(:), intent(in) :: dims
      real, contiguous, dimension(:), pointer, intent(in) :: array_1D
      real, contiguous, dimension(:,:,:,:), pointer, intent(out) :: array_4D

      ! Assign the C pointer to the target Fortran array
      call c_f_pointer(c_loc(array_1D), array_4D, dims)
   end subroutine reshape_real_1D_to_4D

   !> Routine to use c_f_pointer to reshape a contiguous real 1D-array into an 5D-array
   subroutine reshape_real_1D_to_5D(dims, array_1D, array_5D)
      integer, dimension(:), intent(in) :: dims
      real, contiguous, dimension(:), pointer, intent(in) :: array_1D
      real, contiguous, dimension(:,:,:,:,:), pointer, intent(out) :: array_5D

      ! Assign the C pointer to the target Fortran array
      call c_f_pointer(c_loc(array_1D), array_5D, dims)
   end subroutine reshape_real_1D_to_5D

   !> Routine to use c_f_pointer to reshape a contiguous real 1D-array into an 6D-array
   subroutine reshape_real_1D_to_6D(dims, array_1D, array_6D)
      integer, dimension(:), intent(in) :: dims
      real, contiguous, dimension(:), pointer, intent(in) :: array_1D
      real, contiguous, dimension(:,:,:,:,:,:), pointer, intent(out) :: array_6D

      ! Assign the C pointer to the target Fortran array
      call c_f_pointer(c_loc(array_1D), array_6D, dims)
   end subroutine reshape_real_1D_to_6D

   !> Routine to use c_f_pointer to reshape a contiguous real 1D-array into an 7D-array
   subroutine reshape_real_1D_to_7D(dims, array_1D, array_7D)
      integer, dimension(:), intent(in) :: dims
      real, contiguous, dimension(:), pointer, intent(in) :: array_1D
      real, contiguous, dimension(:,:,:,:,:,:,:), pointer, intent(out) :: array_7D

      ! Assign the C pointer to the target Fortran array
      call c_f_pointer(c_loc(array_1D), array_7D, dims)
   end subroutine reshape_real_1D_to_7D

   !> Routine to use c_f_pointer to reshape a contiguous integer 1D-array into a 2D-array
   subroutine reshape_integer_1D_to_2D(dims, array_1D, array_2D)
      integer, dimension(:), intent(in) :: dims
      integer, contiguous, dimension(:), pointer, intent(in) :: array_1D
      integer, contiguous, dimension(:,:), pointer, intent(out) :: array_2D

      ! Assign the C pointer to the target Fortran array
      call c_f_pointer(c_loc(array_1D), array_2D, dims)
   end subroutine reshape_integer_1D_to_2D

   !> Routine to use c_f_pointer to reshape a contiguous integer 1D-array into a 3D-array
   subroutine reshape_integer_1D_to_3D(dims, array_1D, array_3D)
      integer, dimension(:), intent(in) :: dims
      integer, contiguous, dimension(:), pointer, intent(in) :: array_1D
      integer, contiguous, dimension(:,:,:), pointer, intent(out) :: array_3D

      ! Assign the C pointer to the target Fortran array
      call c_f_pointer(c_loc(array_1D), array_3D, dims)
   end subroutine reshape_integer_1D_to_3D

   !> Global index from ndims, dims, indices.
   pure subroutine IDX_XD(ndims, dims, indices, global_index)
      integer, intent(in) :: ndims
      integer, contiguous, dimension(:), intent(in) :: dims, indices
      integer, intent(out) :: global_index

      integer :: d, stride

      global_index = indices(1) - 1
      stride = 1

      do d = 2, ndims
         stride = stride * dims(d - 1)
         global_index = global_index + (indices(d) - 1) * stride
      end do

      ! Adjust index for 1-based indexing of Fortran
      global_index = global_index + 1
   end subroutine IDX_XD

   !> Dimensional indices from ndims, dims, global_index.
   pure subroutine IDX_XD_INV(ndims, dims, global_index, indices)
      integer, intent(in) :: ndims, global_index
      integer, contiguous, dimension(:), intent(in) :: dims
      integer, contiguous, dimension(:), intent(out) :: indices

      integer :: d, remaining_index, stride

      ! Adjust index for 1-based indexing of Fortran
      remaining_index = global_index - 1

      ! Start with the product of all dimensions except the first one
      stride = product(dims(1:ndims-1))

      do d = ndims, 2, -1
         indices(d) = (remaining_index / stride) + 1
         remaining_index = mod(remaining_index, stride)
         stride = stride / dims(d - 1)
      end do

      ! Handle the first dimension separately
      indices(1) = remaining_index + 1

   end subroutine IDX_XD_INV

   !> Routine to print an array of type real
   subroutine print_real_array(ndims, dims, begin_c, end_c, array, stride, title, iounit)
      integer, intent(in) :: ndims, stride, iounit
      integer, dimension(:), intent(in) :: dims, begin_c, end_c
      real, dimension(:), intent(in) :: array
      character(len=*), intent(in) :: title

      integer :: ii, jj, kk, global_index

      character(len=32) :: element_format

      element_format = "(F10.3)"  ! Format for each element

      write(iounit, *) title
      select case(ndims)
       case (1)
         do ii = begin_c(1)+1, end_c(1)
            write(iounit, element_format, advance="no") array(global_index:global_index+stride-1)
         end do

       case (2)
         do ii = begin_c(1)+1, end_c(1)
            write(iounit, *)
            do jj = begin_c(2)+1,end_c(2)
               call IDX_XD(ndims, dims, [ii, jj], global_index)
               write(iounit, element_format, advance="no") array(global_index:global_index+stride-1)
            end do
         end do

       case (3)
         do ii = begin_c(1)+1, end_c(1)
            write(iounit, *)
            write(iounit, *)
            write(iounit, *) "ii: ", ii
            do jj = begin_c(2)+1, end_c(2)
               write(iounit, *)
               do kk = begin_c(3)+1, end_c(3)
                  call IDX_XD(ndims, dims, [ii, jj, kk], global_index)
                  write(iounit, element_format, advance="no") array(global_index:global_index+stride-1)
               end do
            end do
         end do

      end select

   end subroutine print_real_array

   !> Routine to print an array of type integer
   subroutine print_integer_array(ndims, dims, begin_c, end_c, array, stride, title, iounit)
      integer, intent(in) :: ndims, stride, iounit
      integer, dimension(:), intent(in) :: dims, begin_c, end_c
      integer, dimension(:), intent(in) :: array
      character(len=*), intent(in) :: title

      integer :: ii, jj, kk, global_index

      character(len=32) :: element_format

      element_format = "(I6)"  ! Format for each element

      write(iounit, *) title
      select case(ndims)
       case (1)
         do ii = begin_c(1)+1, end_c(1)
            write(iounit, element_format, advance="no") array(global_index:global_index+stride-1)
         end do

       case (2)
         do ii = begin_c(1)+1, end_c(1)
            write(iounit, *)
            do jj = begin_c(2)+1,end_c(2)
               call IDX_XD(ndims, dims, [ii, jj], global_index)
               write(iounit, element_format, advance="no") array(global_index:global_index+stride-1)
            end do
         end do

       case (3)
         do ii = begin_c(1)+1, end_c(1)
            write(iounit, *)
            write(iounit, *)
            write(iounit, *) "ii: ", ii
            do jj = begin_c(2)+1, end_c(2)
               write(iounit, *)
               do kk = begin_c(3)+1, end_c(3)
                  call IDX_XD(ndims, dims, [ii, jj, kk], global_index)
                  write(iounit, element_format, advance="no") array(global_index:global_index+stride-1)
               end do
            end do
         end do

      end select

   end subroutine print_integer_array

   !> Routine to print a 1D array of type integer
   subroutine print_integer_1D_array(array, iounit)
      integer, dimension(:), intent(in) :: array
      integer, intent(in) :: iounit

      integer :: ii

      do ii = 1, size(array)
         write(iounit, "(I6)", advance="no") array(ii)
      end do

   end subroutine print_integer_1D_array

   !> Routine to print a 1D array of type integer
   subroutine print_integer_2D_array(array, iounit)
      integer, dimension(:,:), intent(in) :: array
      integer, intent(in) :: iounit

      integer :: ii, jj

      do ii = 1, size(array, 2)
         write(iounit, *)
         do jj = 1, size(array, 1)
            write(iounit,"(I6)", advance="no") array(jj, ii)
         end do
      end do

   end subroutine print_integer_2D_array

   !> Routine to print a 1D array of type real
   subroutine print_real_1D_array(array, iounit)
      real, dimension(:), intent(in) :: array
      integer, intent(in) :: iounit

      integer :: ii

      do ii = 1, size(array)
         write(iounit, "(F10.3)", advance="no") array(ii)
      end do

   end subroutine print_real_1D_array

   !> Routine to print a 2D array of type real
   subroutine print_real_2D_array(array, iounit)
      real, dimension(:,:), intent(in) :: array
      integer, intent(in) :: iounit

      integer :: ii, jj

      do ii = 1, size(array, 2)
         write(iounit, *)
         do jj = 1, size(array, 1)
            write(iounit,"(F10.3)", advance="no") array(jj, ii)
         end do
      end do

   end subroutine print_real_2D_array

   !> Routine to print a 3D array of type real
   subroutine print_real_3D_array(array, iounit)
      real, dimension(:,:,:), intent(in) :: array
      integer, intent(in) :: iounit

      integer :: ii, jj, kk

      do ii = 1, size(array, 3)
         write(iounit, *)
         write(iounit, *)
         write(iounit, *) "ii: ", ii
         do jj = 1, size(array, 2)
            write(iounit, *)
            do kk = 1, size(array, 1)
               write(iounit,"(F10.3)", advance="no") array(kk, jj, ii)
            end do
         end do
      end do

   end subroutine print_real_3D_array

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

      dx = abs((domain_end - domain_begin) / (grid_size_with_ghosts - 1))
   end subroutine calculate_dx

   !> Routine to calculate Courant-Friedrichs-Lewy (CFL) number
   pure subroutine calculate_CFL(magnitude_velocity, dx, dt, CFL)
      real, dimension(:), intent(in) :: magnitude_velocity, dx
      real, intent(in) :: dt
      real, intent(out) :: CFL

      CFL = dt * sum(magnitude_velocity / dx)
   end subroutine calculate_CFL

   !> Routine to find dt from CFL number
   pure subroutine calculate_dt_from_CFL(CFL, magnitude_velocity, dx, dt)
      real, intent(in) :: CFL
      real, dimension(:), intent(in) :: magnitude_velocity, dx
      real, intent(out) :: dt

      dt = CFL / sum(magnitude_velocity / dx)
   end subroutine calculate_dt_from_CFL

   !> Subroutine to swap two pointers
   pure subroutine swap_pointers(ptr1, ptr2)
      real, contiguous, dimension(:), pointer, intent(inout) :: ptr1, ptr2
      real, contiguous, dimension(:), pointer :: temp_ptr

      temp_ptr => ptr1
      ptr1 => ptr2
      ptr2 => temp_ptr
   end subroutine swap_pointers

   !> Subroutine to swap two 2D-pointers
   pure subroutine swap_pointers_2D(ptr1, ptr2)
      real, contiguous, dimension(:,:), pointer, intent(inout) :: ptr1, ptr2
      real, contiguous, dimension(:,:), pointer :: temp_ptr

      temp_ptr => ptr1
      ptr1 => ptr2
      ptr2 => temp_ptr
   end subroutine swap_pointers_2D

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
