      module mod_integration

      use mod_common
      use mod_structure
      use mod_recherche
      use mod_diffusion
      use mod_reaction
      use mod_arbre

      contains

      include "integration_reaction.f90"
      include "integration_diffusion.f90"
      include "construction_phi.f90"

      end module mod_integration
