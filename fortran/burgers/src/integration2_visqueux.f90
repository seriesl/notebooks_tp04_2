      subroutine integration2_visqueux(dir)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: dir

      INTEGER :: l, l_dt, leaf, m
      INTEGER, DIMENSION(4,3) :: ijk_dt

      REAL(kind=8) :: sigma, rho_ec

      REAL(kind=8), DIMENSION(nvar) :: fvisc_gauche, fvisc_droit

      TYPE(structure_maille), pointer :: maille

      TYPE(structure_maille), pointer ::  local

! .....definition des indices des cellules a droite :
!      ijk_dt( 2**(ndim-1) , dir )
  
      ijk_dt(1,1) = 2
      ijk_dt(2,1) = 4
      ijk_dt(3,1) = 6
      ijk_dt(4,1) = 8
  
      ijk_dt(1,2) = 3
      ijk_dt(2,2) = 4
      ijk_dt(3,2) = 7
      ijk_dt(4,2) = 8
  
      ijk_dt(1,3) = 5
      ijk_dt(2,3) = 6
      ijk_dt(3,3) = 7
      ijk_dt(4,3) = 8

! *****On parcourt les feuilles reelles (non-fictives)

      Do leaf = 1, compteur_feuille

        maille => feuille_reelle(leaf)%ptr

! =================== Integration Predicteur  Visqueux =============

! .....Rapport : sigma = dt / dx
        sigma = dt / maille%dx(dir)

! +++++Maille courante

        local => maille

! =====Flux visqueux a droite

        fvisc_droit(:) = maille%fvisc_dt(:,dir)

! =====Flux visqueux a gauche

! .....par defaut, on garde le flux a gauche 
        fvisc_gauche(:) = maille%fvisc_gch(:,dir)

! .....Si la feuille precedente existe
        if( associated(maille%support(dir)%precedent) ) then

          local => maille%support(dir)%precedent

! .....si la feuille precedente est reelle
          if( local%feuille ) then

! .....on transfert le flux a droite pour conservation
            fvisc_gauche(:) = local%fvisc_dt(:,dir)

          endif

! .....Si la feuille precedente est a un niveau plus eleve

          if( associated(local%fils(1)%ptr) ) then

! .....si les fils de la feuille precedente sont reels
            if( local%fils(1)%ptr%feuille ) then

! .....Somme des flux Visqueux pour conservation

              fvisc_gauche(:) = 0.

              do l = 1, 2**(ndim-1)

                fvisc_gauche(:) = fvisc_gauche(:) + &
                          local%fils( ijk_dt(l,dir) )%ptr%fvisc_dt(:,dir)

              enddo

              fvisc_gauche(:) = fvisc_gauche(:) / float(2**(ndim-1))

            endif

          endif

        else
   
          if( maille%maille_i(dir) /= 1 ) then
            Print*, ' Integration - Pb de CAL < : maille niv= ', &
             maille%maille_niv, ' i = ', maille%maille_i
            Stop
          
          endif

        endif

! *****Predicteur Flux visqueux

        do m = 1, nvar
          maille%unp1(m) = 0.5*( maille%u(m) + maille%u0(m) ) - &
                      0.5*sigma*( fvisc_droit(m) - fvisc_gauche(m) )
        enddo

      Enddo

      end subroutine integration2_visqueux
