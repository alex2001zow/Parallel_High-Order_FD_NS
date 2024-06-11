module constants_module
   implicit none

   private

   integer, parameter ::character_len = 255

   integer, parameter :: MASTER_RANK = 0

   integer, parameter :: neighbor_current_rank = -10 ! Replace with MPI_PROC_NULL
   integer, parameter :: neighbor_non_existant_rank = -1 ! Fix spelling...

   real, parameter :: pi = 3.1415926535897932384626433832795

   public :: MASTER_RANK, neighbor_current_rank, neighbor_non_existant_rank
   public :: pi

end module constants_module
