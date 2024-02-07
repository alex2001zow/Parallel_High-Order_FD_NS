module constants_module
   implicit none

   private

   integer, parameter :: MASTER_RANK = 0

   integer, parameter :: neighbor_current_rank = -10
   integer, parameter :: neighbor_non_existant_rank = -1

   public :: MASTER_RANK, neighbor_current_rank, neighbor_non_existant_rank

end module constants_module