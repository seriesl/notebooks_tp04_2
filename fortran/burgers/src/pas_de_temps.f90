	subroutine pas_de_temps

	use mod_common
	use mod_structure
        use mod_fonction

	implicit none

        INTEGER :: l, leaf, m

        REAL(kind=8) :: vp, mu_0

	TYPE(structure_maille), pointer :: maille

! +++++On parcourt les feuilles reelles (non-fictives)

        Do leaf = 1, compteur_feuille

          maille => feuille_reelle(leaf)%ptr

! .....Rayon spectral
          vp = 0.
          do m = 1,nvar
            vp = max( vp, abs(maille%u(m)) )
          enddo

! .....Vicosite moleculaire et conductivite thermique
          call viscosite(maille, mu_0)

! .....Calcule du pas de temps avec critere CFL
          if( Reynolds /= 0. ) then
            do l =1, ndim
              dt = min( dt, maille%dx(l)*CFL/vp, &
                    maille%dx(l)**2*Reynolds*.90*CFL / mu_0 )
            enddo
          else
            do l =1, ndim
              dt = min( dt, maille%dx(l)*CFL/vp )
            enddo
          endif

	Enddo

	end subroutine pas_de_temps
