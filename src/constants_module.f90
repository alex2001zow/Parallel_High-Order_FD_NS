module constants_module
   implicit none

   private

   integer, parameter ::character_len = 255

   character(character_len), parameter :: filename_txt = "output/output_from_" // repeat(" ", &
      character_len - len_trim("output/output_from_"))
   character(character_len), parameter :: filename_dat = "output/system_solution.dat" // repeat(" ", &
      character_len - len_trim("output/system_solution.dat"))
   character(character_len), parameter :: filename_layout_dat = "output/system_solution_layout.dat" // repeat(" ", &
      character_len - len_trim("output/system_solution_layout.dat"))

   integer, parameter :: MASTER_RANK = 0

   integer, parameter :: neighbor_current_rank = -10
   integer, parameter :: neighbor_non_existant_rank = -1

   real, parameter :: pi = 3.1415926535897932384626433832795

   public :: filename_txt, filename_dat, filename_layout_dat
   public :: MASTER_RANK, neighbor_current_rank, neighbor_non_existant_rank
   public :: pi

end module constants_module
