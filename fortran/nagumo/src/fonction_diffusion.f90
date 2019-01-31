      subroutine fonction_diffusion(ni, phi, dphidt, rpar, ipar)

      use mod_common
      use mod_structure
      use mod_arbre

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ni, ipar ! -> numero de variable

      DOUBLE PRECISION,INTENT(INOUT),DIMENSION(ni) :: phi 

      DOUBLE PRECISION,INTENT(OUT),DIMENSION(ni) :: dphidt
 
      DOUBLE PRECISION, INTENT(IN) :: rpar


      INTEGER :: ir, jr, kr, comp_f
      INTEGER :: inc, k, l, m, n, dir, i

      TYPE(structure_maille), POINTER :: courant

! .....Mise a jour des valeurs dans les cellules

      comp_f = 0
     
       do i = 1,nvar ! supprimer si var separees
          do l = 1,compteur_feuille
!!             tab_phi(l)%var(ipar)%t_ptr = phi(l) ! var separees
             tab_phi(l)%var(i)%t_ptr = phi((i-1)*compteur_feuille + l)
          enddo
       enddo ! supprimer si var separees


 
! .....Codage par valeurs moyennes 

              do kr = 1, nbre_racine(3)
                do jr = 1, nbre_racine(2)
                  do ir = 1, nbre_racine(1)
  
                    courant => racine(ir, jr, kr)%ptr
                    call encoding_moy_diff(courant,ipar) !var separees

                  enddo
                enddo
             enddo

! .....Valorisation des feuilles fictives

                    call valeur_feuille_fictive(ipar) !var separees
                    
! .....Calcul des flux diffusifs aux interfaces

                    call flux_diff_1(ipar) ! var separees

! .....Calcul de dphidt

                   call flux_diff_2(ni, dphidt, ipar) ! var separees
   
      end subroutine fonction_diffusion
