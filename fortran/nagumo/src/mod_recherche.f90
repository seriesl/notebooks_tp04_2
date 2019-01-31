      module mod_recherche

      use mod_common
      use mod_structure

      contains

      include "recherche_support.f90"
      include "recherche_support_diff.f90"
      include "cal_indice.f90"
      include "creation.f90"
      include "recherche.f90"

      end module mod_recherche
