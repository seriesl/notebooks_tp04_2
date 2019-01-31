!	module de l'integration du terme reactif (RADAU5)

      module mod_reaction

      use mod_common

      implicit none

      INTEGER  :: ITOLREAC
      INTEGER  :: IJOMREAC
      INTEGER  :: MLJACREAC 
      INTEGER  :: MUJACREAC 
      INTEGER  :: IMASREAC
      INTEGER  :: MLMASREAC
      INTEGER  :: MUMASREAC
      INTEGER  :: IPARREAC
      INTEGER  :: IOUTREAC
      INTEGER  :: LENWKREAC
      INTEGER  :: LENIWKREAC

      DOUBLE PRECISION :: HRREAC
      DOUBLE PRECISION, dimension(ndim) :: RPARREAC 
      DOUBLE PRECISION, dimension(:), allocatable :: ELWRKREAC
      
      INTEGER, dimension(:), allocatable :: IELWRKREAC

      INTEGER :: ISTATEREAC


      contains


  
      include "reaceq.f90"
 
      end module mod_reaction
