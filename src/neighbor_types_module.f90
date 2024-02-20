module neighbor_types_module
   implicit none

   private

   !> Enumerators for the different types of neighbors in 1D
   enum, bind(c)
      enumerator :: LEFT_1D = 1, CENTER_1D = 2, RIGHT_1D = 3
   end enum

   !> Enumerators for the different types of neighbors in 2D
   enum, bind(c)
      enumerator :: TOP_LEFT_2D = 1, TOP_2D = 2, TOP_RIGHT_2D = 3
      enumerator :: LEFT_2D = 4, CENTER_2D = 5, RIGHT_2D = 6
      enumerator :: BOTTOM_LEFT_2D = 7, BOTTOM_2D = 8, BOTTOM_RIGHT_2D = 9
   end enum

   public ::get_neighbor_range

contains

   !> Subroutine to get the range per dimension for a given neighbor in 1D, 2D and one day 3D. DOUBLE CHECK THIS ROUTINE!
   subroutine get_neighbor_range(ndims, neighbor, dims, begin, end)
      integer, intent(in) :: ndims, neighbor
      integer, dimension(ndims), intent(in) :: dims
      integer, dimension(ndims), intent(out) :: begin, end

      if (ndims == 1) then
         select case(neighbor)
          case(LEFT_1D)
            begin(1) = 1
            end(1) = 1
          case(CENTER_1D)
            begin(1) = 2
            end(1) = dims(1) - 1
          case(RIGHT_1D)
            begin(1) = dims(1)
            end(1) = dims(1)
          case default
            print *, "Invalid neighbor type"
            begin(1) = -1
            end(1) = -1
         end select
      end if

      if (ndims == 2) then
         select case(neighbor)
          case(TOP_LEFT_2D)
            begin = [1, 1]
            end = [1, 1]
          case(TOP_2D)
            begin = [1, 2]
            end = [1, dims(2) - 1]
          case(TOP_RIGHT_2D)
            begin = [1, dims(2)]
            end = [1, dims(2)]
          case(LEFT_2D)
            begin = [2, 1]
            end = [dims(1) - 1, 1]
          case(CENTER_2D)
            begin = [2, 2]
            end = [dims(1) - 1, dims(2) - 1]
          case(RIGHT_2D)
            begin = [2, dims(2)]
            end = [dims(1) - 1, dims(2)]
          case(BOTTOM_LEFT_2D)
            begin = [dims(1), 1]
            end = [dims(1), 1]
          case(BOTTOM_2D)
            begin = [dims(1), 2]
            end = [dims(1), dims(2) - 1]
          case(BOTTOM_RIGHT_2D)
            begin = [dims(1), dims(2)]
            end = [dims(1), dims(2)]
          case default
            print *, "Invalid neighbor type"
            begin = [-1, -1]
            end = [-1, -1]
         end select
      end if
   end subroutine get_neighbor_range


end module neighbor_types_module
