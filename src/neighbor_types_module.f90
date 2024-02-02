module neighbor_types_module
   implicit none

   private

   !> Type for storing the range of indices for a given neighbor
   type index_range
      integer :: ndims
      integer, allocatable :: begin(:), end(:)
   end type index_range

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

   public :: LEFT_1D, CENTER_1D, RIGHT_1D, TOP_LEFT_2D, TOP_2D, TOP_RIGHT_2D, LEFT_2D, CENTER_2D, RIGHT_2D, BOTTOM_LEFT_2D, BOTTOM_2D, BOTTOM_RIGHT_2D

end module neighbor_types_module
