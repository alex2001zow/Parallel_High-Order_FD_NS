module constants_module
   implicit none

   private

   ! Set output filename
   character(255), parameter :: filename = "output/output_from_" // repeat(" ", 255 - len_trim("output/output_from_"))

   integer, parameter :: MASTER_RANK = 0

   integer, parameter :: neighbor_current_rank = -10
   integer, parameter :: neighbor_non_existant_rank = -1

   real, parameter :: pi = 3.1415926535897932384626433832795

   public :: filename, MASTER_RANK, neighbor_current_rank, neighbor_non_existant_rank, pi

end module constants_module
