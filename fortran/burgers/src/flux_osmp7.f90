      subroutine flux_osmp7(dir)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: dir

      LOGICAL :: flag
 
      INTEGER :: i_sym, l, l_dt, leaf, m
      INTEGER :: s_deb, s_fin
      INTEGER :: ijk_gch
      INTEGER, DIMENSION(4,3) :: ijk_dt

      REAL(kind=8), DIMENSION(nvar, -nordre+1:nordre) :: u_n
      REAL(kind=8), DIMENSION(nvar) :: C_mp
      REAL(kind=8), DIMENSION(nvar) :: fmp_droit

      REAL(kind=8), DIMENSION(ndim) :: x_0

      REAL(kind=8) :: rayon
      REAL(kind=8) :: rho_ec, sigma
      REAL(kind=8) :: rho_wall, rhoE_wall

      TYPE(structure_maille), pointer :: maille
      TYPE(structure_maille), pointer :: local, courant

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

! ======== Correction OSMP7 : Condition a la limite amont du domaine =====

        if(maille%maille_i(dir) == 1 .and. caltype_deb(dir) /= 2) then

          local => maille

! .....Rapport : sigma = dt / dx
          sigma = dt / local%dx(dir)

! +++++Pour toutes les mailles utiles

          s_deb = 1
          s_fin = 1

          do l = s_deb, nordre

            do m = 1, nvar
              u_n(m,l) = local%u0(m)
            enddo

! .....on passe a la maille suivante
            s_fin = l

!!CT0110            if( .not.associated(local%support(dir)%suivant) ) exit
            local => local%support(dir)%suivant

          enddo

! +++++Reconstruction en amont du domaine 

! .....Conditions de symetrie

          if( caltype_deb(dir) == -1 ) then

            do l = -nordre+1, s_deb-1

              i_sym = 2*s_deb-(l+1)

              do m = 1, nvar
                u_n(m,l) = u_n(m,i_sym)
              enddo
              u_n(dir+1,l) = -u_n(dir+1,i_sym)

            enddo

! .....Conditions de frontiere fluide ==> derivee nulle

          elseif( caltype_deb(dir) == 0 ) then

            do l = -nordre+1, s_deb-1

              i_sym = 2*s_deb-(l+1)

              do m = 1, nvar
!!                u_n(m,l) = u_n(m,i_sym)
                u_n(m,l) = u_n(m,s_deb)
              enddo

            enddo

! .....Conditions de paroi

          elseif( caltype_deb(dir) == 1 ) then

            do l = -nordre+1, s_deb-1
              do m = 1, nvar
                u_n(m, l) = 0.
              enddo
            enddo
            
          endif

! +++++Calcul de la correction OSMP7 a la frontiere amont du domaine

          call osmp7(dir, sigma, u_n, C_mp)

! .....Correction TVD-MP a gauche de la maille courante

          maille%fmp_gch(:) = C_mp(:)

        endif

! ============ Correction OSMP7 a droite pour toutes les feuilles ========

        flag = .false.

! .....Test pour reconstruire niveau feuille ou niveau sup ?

        if( associated(maille%support(dir)%suivant) ) then

          if( associated(maille%support(dir)%suivant%fils(1)%ptr) ) then

            if(maille%support(dir)%suivant%fils(1)%ptr%feuille) &
              flag = .true.

          endif
        endif

! *****Si la feuille a des fils et que 
!      la feuille suivante a des fils non-fictifs

        if( flag .and. associated(maille%fils(1)%ptr) ) then

          fmp_droit(:) = 0.
 
! +++++La correction est evaluee au niveau superieur (Fils de la feuille)

          do l_dt = 1, 2**(ndim-1)

            local => maille%fils( ijk_dt(l_dt,dir) )%ptr
            courant => maille%fils( ijk_dt(l_dt,dir) )%ptr

! .....Rapport : sigma = dt / dx
            sigma = dt / local%dx(dir)

! *****initialisation pour toutes les mailles utiles

            do l = 0, -nordre+1, -1
              s_deb = l
              if( .not.associated(local%support(dir)%precedent) &
                  .or. l == -nordre+1 ) exit
              local => local%support(dir)%precedent
            enddo

! .....Pour toutes les mailles utiles

            s_fin = s_deb
            do l = s_deb, nordre

              do m = 1, nvar
                u_n(m,l) = local%u0(m)
              enddo

! .....on passe a la maille suivante
              s_fin = l

              if( .not.associated(local%support(dir)%suivant) ) exit
              local => local%support(dir)%suivant

            enddo

! *****Conditions a la limite amont

            if( s_deb > -nordre+1 ) then

! .....Conditions de symetrie

              if( caltype_deb(dir) == -1 ) then

                do l = -nordre+1, s_deb-1

                  i_sym = 2*s_deb-(l+1)
 
                  do m = 1, nvar
                    u_n(m,l) = u_n(m,i_sym)
                  enddo
                  u_n(dir+1,l) = -u_n(dir+1,i_sym)

                enddo

! .....Conditions de frontiere fluide ==> derivee nulle

              elseif( caltype_deb(dir) == 0 ) then

                do l = -nordre+1, s_deb-1

                  i_sym = 2*s_deb-(l+1)

                  do m = 1, nvar
!!                    u_n(m,l) = u_n(m,i_sym)
                    u_n(m,l) = u_n(m,s_deb)
                  enddo

                enddo

! .....Conditions de paroi

              elseif( caltype_deb(dir) == 1 ) then

                do l = -nordre+1, s_deb-1
                  do m = 1, nvar
                    u_n(m, l) = 0.
                  enddo
                enddo
            
              endif

            endif

! *****Conditions a la limite aval

            if( s_fin < nordre ) then

! .....Conditions de symetrie

              if( caltype_fin(dir) == -1 ) then

                do l = s_fin+1, nordre

                  i_sym = 2*s_fin-(l-1)

                  do m = 1, nvar
                    u_n(m,l) = u_n(m,i_sym)
                  enddo
                  u_n(dir+1,l) = -u_n(dir+1,i_sym)

                enddo

! .....Conditions de frontiere fluide ==> derivee nulle

              elseif( caltype_fin(dir) == 0 ) then

                do l = s_fin+1, nordre

                  i_sym = 2*s_fin-(l-1)

                  do m = 1, nvar
!!                    u_n(m,l) = u_n(m,i_sym)
                    u_n(m,l) = u_n(m,s_fin)
                  enddo

                enddo

! .....Conditions de paroi

              elseif( caltype_fin(dir) == 1 ) then

                do l = s_fin+1, nordre
                  do m = 1, nvar
                    u_n(m, l) = 0.
                  enddo
                enddo
            
              endif

            endif

! +++++Calcul de la correction OSMP7 a droite

            call osmp7(dir, sigma, u_n, C_mp)

! .....Somme des flux pour conservation

            fmp_droit(:) = fmp_droit(:) + &
                            C_mp(:) / float(2**(ndim-1))

! +++++Transfert des flux pour le flux a gauche de la feuille suivante 

! .....indice de la maille suivante
            ijk_gch = ijk_dt( l_dt,dir ) - 2**(dir-1)

            maille%support(dir)%suivant%fils(ijk_gch)%ptr%fmp_gch(:)=C_mp(:)

          enddo

! .....Correction TVD-MP a droite pour la maille courante

          maille%fmp_dt(:) = fmp_droit(:) 

        else

! +++++Sinon, on reconstruit au niveau de la feuille

          local => maille
          courant => maille

! .....Rapport : sigma = dt / dx
          sigma = dt / local%dx(dir)

! .....initialisation du support

          do l = 0, -nordre+1, -1
            s_deb = l
            if( .not.associated(local%support(dir)%precedent) &
                .or. l == -nordre+1 ) exit
            local => local%support(dir)%precedent
          enddo

! *****Pour toutes les mailles utiles

          s_fin = s_deb
          do l = s_deb, nordre

            do m = 1, nvar
              u_n(m,l) = local%u0(m)
            enddo

! .....on passe a la maille suivante
            s_fin = l

            if( .not.associated(local%support(dir)%suivant) ) exit
            local => local%support(dir)%suivant

          enddo

! *****Conditions a la limite amont

          if( s_deb > -nordre+1 ) then

! .....Conditions de symetrie

            if( caltype_deb(dir) == -1 ) then

              do l = -nordre+1, s_deb-1

                i_sym = 2*s_deb-(l+1)

                do m = 1, nvar
                  u_n(m,l) = u_n(m,i_sym)
                enddo
                u_n(dir+1,l) = -u_n(dir+1,i_sym)

              enddo

! .....Conditions de frontiere fluide ==> derivee nulle

            elseif( caltype_deb(dir) == 0 ) then

              do l = -nordre+1, s_deb-1

                i_sym = 2*s_deb-(l+1)

                do m = 1, nvar
!!                  u_n(m,l) = u_n(m,i_sym)
                  u_n(m,l) = u_n(m, s_deb)
                enddo

              enddo

! .....Conditions de paroi

            elseif( caltype_deb(dir) == 1 ) then

              do l = -nordre+1, s_deb-1
                do m = 1, nvar
                  u_n(m, l) = 0.
                enddo
              enddo
            
            endif

          endif

! *****Conditions a la limite aval

          if( s_fin < nordre ) then

! .....Conditions de symetrie

            if( caltype_fin(dir) == -1 ) then

              do l = s_fin+1, nordre

                i_sym = 2*s_fin-(l-1)

                do m = 1, nvar
                  u_n(m,l) = u_n(m,i_sym)
                enddo
                u_n(dir+1,l) = -u_n(dir+1,i_sym)

              enddo

! .....Conditions de frontiere fluide ==> derivee nulle

            elseif( caltype_fin(dir) == 0 ) then

              do l = s_fin+1, nordre

                i_sym = 2*s_fin-(l-1)

                do m = 1, nvar
!!                  u_n(m,l) = u_n(m,i_sym)
                  u_n(m,l) = u_n(m, s_fin)
                enddo

              enddo

! .....Conditions de paroi

            elseif( caltype_fin(dir) == 1 ) then

              do l = s_fin+1, nordre
                do m = 1, nvar
                  u_n(m, l) = 0.
                enddo
              enddo
            
            endif

          endif

! +++++Calcul de la correction OSMP7 a droite

          call osmp7(dir, sigma, u_n, C_mp)

! .....Correction TVD-MP a droite pour la maille courante

          maille%fmp_dt(:) = C_mp(:)

! +++++Transfert des flux pour le flux a gauche de la feuille suivante 

          if( associated(maille%support(dir)%suivant) ) &
            maille%support(dir)%suivant%fmp_gch(:) = C_mp(:) 

        endif

      Enddo

      end subroutine flux_osmp7
