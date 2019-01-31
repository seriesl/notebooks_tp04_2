      subroutine flux_diff_1(mm)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      INTEGER :: dir

      LOGICAL :: flag
 
      INTEGER :: i_dir, i_sym, l, l_dt, m, mm, ni, feu
      INTEGER :: s_deb, s_fin
      INTEGER :: ijk_gch
      INTEGER, DIMENSION(3) :: i_beg, i_end, i_pas
      INTEGER, DIMENSION(4,3) :: ijk_dt

      DOUBLE PRECISION, DIMENSION(nvar) :: u_0
      DOUBLE PRECISION, DIMENSION(nvar) :: u_1


      DOUBLE PRECISION, DIMENSION(nvar) :: fdiff_droit, C_diff

      TYPE(structure_maille), pointer :: courant, crt_suivant, local


! .....Definition des indices des cellules a droite : 
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


!..... On parcours le tableau de feuilles

      do feu = 1, compteur_feuille

     
! .....Test pour reconstruire au niveau de la feuille ou au niveau superieur ?


      do dir = 1, ndim

        flag = .false.

        if( associated(tab_phi(feu)%support(dir)%suivant) ) then

            if( associated(tab_phi(feu)%support(dir)%suivant%fils(1)%ptr) ) then

            if(tab_phi(feu)%support(dir)%suivant%fils(1)%ptr%feuille) &
              flag = .true.

          endif
        endif

! *****Si la feuille a des fils (qui seront fictives) et que 
!      la feuille suivante a des fils non-fictifs

        if( flag .and. associated(tab_phi(feu)%fils(1)%ptr) ) then

! ======== Flux visqueux a droite pour toutes les feuilles =====

          fdiff_droit(:) = 0.d0

! +++++La correction est evaluee au niveau superieur (Fils de la feuille)

          do l_dt = 1, 2**(ndim-1)

! *****Maille courante

            courant => tab_phi(feu)%fils( ijk_dt(l_dt,dir) )%ptr

! .....vitesses au centre de la maille
              do m = 1, nvar
                u_0(m) = courant%u(m)
              enddo

! *****Maille suivante

            if( associated(courant%support(dir)%suivant) ) then

              crt_suivant => courant%support(dir)%suivant

! .....vitesses au centre de la maille
                do m = 1, nvar
                  u_1(m) = crt_suivant%u(m)
                enddo
! ++++++Calcule des flux

            C_diff(:) = ( u_1(:) - u_0(:) ) / courant%dx(dir)


! +++++Somme des flux pour conservation

            fdiff_droit(:) = fdiff_droit(:) + C_diff(:) / dble(2**(ndim-1))

! +++++Transfert des flux pour le flux a gauche de la feuille suivante 

                crt_suivant%fdiff_gch(:,dir) = C_diff(:)


            else
                do m = 1, nvar
                  u_1(m) = u_0(m)
                enddo


! ++++++Calcule des flux

            C_diff(:) = ( u_1(:) - u_0(:) ) / courant%dx(dir)


! +++++Somme des flux pour conservation

            fdiff_droit(:) = fdiff_droit(:) + C_diff(:) / dble(2**(ndim-1))

            endif

          enddo

! +++++Sinon, on reconstruit au niveau de la feuille

        else

! ======== Flux visqueux a droite pour toutes les feuilles =====

          fdiff_droit(:) = 0.d0


! .....valeurs au centre de la maille

            do m = 1, nvar
              u_0(m) = tab_phi(feu)%var(m)%t_ptr
            enddo

! *****Maille suivante

          if( associated(tab_phi(feu)%support(dir)%suivant) ) then

            crt_suivant => tab_phi(feu)%support(dir)%suivant


! .....valeurs au centre de la maille
              do m = 1, nvar
                u_1(m) = crt_suivant%u(m)
              enddo
! ++++++Calcule des flux

            C_diff(:) = ( u_1(:) - u_0(:) ) / tab_phi(feu)%dx(dir)


! +++++"Somme"des flux pour conservation

          fdiff_droit(:) = C_diff(:)

! +++++Transfert du flux pour le flux a gauche de la feuille suivante,
!      si elle existe 

            crt_suivant%fdiff_gch(:,dir) = C_diff(:) 


          else
              do m = 1, nvar
                u_1(m) = u_0(m)
              enddo

! ++++++Calcule des flux

            C_diff(:) = ( u_1(:) - u_0(:) ) / tab_phi(feu)%dx(dir)


! +++++"Somme"des flux pour conservation

          fdiff_droit(:) = C_diff(:)

         endif

        endif

! .....Flux Visqueux a droite de la maille courante

        do m = 1,nvar
           tab_phi(feu)%fdiff_dt(m,dir)%t_ptr = fdiff_droit(m)
        enddo 

      enddo 



      enddo

      end subroutine flux_diff_1
