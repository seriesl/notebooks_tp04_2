      subroutine flux_diff_2(ni,dphidt,mm)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER :: dir, ni, mm

      DOUBLE PRECISION,INTENT(INOUT),DIMENSION(ni) :: dphidt 

      INTEGER :: l, l_dt, m, feu
      INTEGER, DIMENSION(4,3) :: ijk_dt

      DOUBLE PRECISION :: sigma, rho_ec

      DOUBLE PRECISION, DIMENSION(nvar) :: fdiff_gauche, fdiff_droit

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


!..... On parcours le tableau de feuilles

      do feu = 1, compteur_feuille

       do m = 1, nvar

          dphidt( (m-1)*compteur_feuille + feu ) = 0.d0

       enddo

      do dir = 1, ndim

! =================== Integration Predicteur  Visqueux =============

! .....Rapport : sigma = 1 / dx
        sigma = 1. / tab_phi(feu)%dx(dir)

! +++++Maille courante

! =====Flux visqueux a droite

        do m = 1,nvar
           fdiff_droit(m) = tab_phi(feu)%fdiff_dt(m,dir)%t_ptr
        enddo


! =====Flux visqueux a gauche

! .....par defaut, on garde le flux a gauche 
        do m = 1,nvar
           fdiff_gauche(m) = tab_phi(feu)%fdiff_gch(m,dir)%t_ptr
        enddo

! .....Si la feuille precedente existe
        if( associated(tab_phi(feu)%support(dir)%precedent) ) then

          local => tab_phi(feu)%support(dir)%precedent

! .....Si la feuille precedente est a un niveau plus eleve

          if( associated(local%fils(1)%ptr) ) then

! .....si les fils de la feuille precedente sont reels
            if( local%fils(1)%ptr%feuille ) then

! .....Somme des flux Visqueux pour conservation

              fdiff_gauche(:) = 0.

              do l = 1, 2**(ndim-1)

                fdiff_gauche(:) = fdiff_gauche(:) + &
                local%fils( ijk_dt(l,dir) )%ptr%fdiff_dt(:,dir)/dble(2**(ndim-1))

              enddo


            endif

          endif

        else
 
! *****Conditions a la limite amont dans la direction dir


              fdiff_gauche(:) = 0.d0

        endif


! *****Calcule Dd2f/dx2

        do m = 1, nvar

          dphidt( (m-1)*compteur_feuille + feu ) = &
          dphidt( (m-1)*compteur_feuille + feu ) + &
          Dif(m)*sigma*(fdiff_droit(m) - fdiff_gauche(m))

 
       enddo

      enddo 


      enddo


      end subroutine flux_diff_2
