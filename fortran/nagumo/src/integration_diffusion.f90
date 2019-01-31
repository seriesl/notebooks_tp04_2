      subroutine integration_diffusion

      use mod_common
      use mod_structure
      use mod_diffusion

      IMPLICIT NONE

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: p_phi

      DOUBLE PRECISION :: x, xend, rpardiff
      DOUBLE PRECISION, dimension(:), allocatable :: ELWRKDIFF
      INTEGER :: n, l, i, ipardiff

       n = nvar*compteur_feuille
!!       n = compteur_feuille ! var separees
      
       allocate (p_phi(n))

       allocate (ELWRKDIFF(8*n))

!!      DO i = 1,nvar ! var separees
 
       do i = 1,nvar ! supprimer si var separees
          do l = 1,compteur_feuille
!!             p_phi(l) = tab_phi(l)%var(i)%t_ptr ! var separees
             p_phi((i-1)*compteur_feuille + l) = tab_phi(l)%var(i)%t_ptr
          enddo
       enddo ! supprimer si var separees


       rpardiff = 0.
       ipardiff = i ! pour var separees

       x = temps
       xend = x + dt
       atoldiff = 1.d-10
       rtoldiff = 1.d-10
       !!atoldiff = 1.d-15
       !!rtoldiff = 1.d-15
       !!HRDIFF = dt
       HRDIFF = 1.e-6

      CALL ROCK4(n,x,xend,HRDIFF,p_phi,diffeq,ATOLDIFF,&
                 RTOLDIFF,ELWRKDIFF,IELWRKDIFF,ISTATEDIFF,&
                 RPARDIFF,IPARDIFF)



       do i = 1,nvar ! supprimer si var separees
          do l = 1,compteur_feuille
             if (p_phi((i-1)*compteur_feuille + l) < 0D0) &
                 p_phi((i-1)*compteur_feuille + l) = 0d0
!!             tab_phi(l)%var(i)%t_ptr = p_phi(l) ! var separees
             tab_phi(l)%var(i)%t_ptr = p_phi((i-1)*compteur_feuille + l)
          enddo
       enddo ! supprimer si var separees

!!      ENDDO ! var separees


      deallocate(ELWRKDIFF,p_phi)


      end subroutine integration_diffusion
