      subroutine residu

      use mod_common
      use mod_structure

      implicit none

      INTEGER :: l, leaf, m

      REAL(kind=8) :: dvol
      REAL(kind=8), DIMENSION(nvar) :: du

      TYPE(structure_maille), pointer :: maille

! *****On parcourt les feuilles reelles (non-fictives)

      Do leaf = 1, compteur_feuille

        maille => feuille_reelle(leaf)%ptr

! +++++calcul des residus

        do m = 1, nvar
          du(m) = maille%u(m) - maille%u_ndt(m)
        enddo

! +++++Volume de la maille

        dvol = 1.
        do l = 1, ndim
          dvol = dvol * maille%dx(l)
        enddo

! +++++Volume total
        Volume = Volume + dvol

! +++++calcul norme L2 des residus
 
        do m = 1, nvar

          if( abs(du(m)) >= err_max(m) ) then
            err_max(m) = abs(du(m))
            ijk_resl2max(:, m) = maille%maille_i(:)
          endif

!!          errol2(m) = errol2(m) + (du(m)/dt)**2
          errol2(m) = errol2(m) + du(m)**2

        enddo

! *****evaluation des integrales pour interaction spot-choc

        integral(:) = integral(:) + maille%u(:) * dvol

      Enddo

      end subroutine residu
