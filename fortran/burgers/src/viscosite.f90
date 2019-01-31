      subroutine viscosite(maille, mu_loc)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      REAL(kind=8), INTENT(OUT) :: mu_loc

      INTEGER :: l, m

      REAL(kind=8) :: ec, T_loc

! +++++Vicosite dynamique

      mu_loc = 1.

      end subroutine viscosite
