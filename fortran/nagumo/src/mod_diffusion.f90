!	module de l'integration du terme diffusif (ROCK4)

      module mod_diffusion

      implicit none

      INTEGER, PARAMETER :: LENIWKDIFF = 12

      DOUBLE PRECISION :: HRDIFF 
      DOUBLE PRECISION :: ATOLDIFF 
      DOUBLE PRECISION :: RTOLDIFF 

    
      INTEGER, dimension(LENIWKDIFF) :: IELWRKDIFF

      INTEGER :: ISTATEDIFF

      contains

      include "diffeq.f90"
      include "flux_diff_1.f90"
      include "flux_diff_2.f90"
      include "fonction_diffusion.f90"

      end module mod_diffusion
