      subroutine integration_reaction

      use mod_common
      use mod_structure
      use mod_reaction

      IMPLICIT NONE

 
      DOUBLE PRECISION :: x, xend, atol, rtol 
      DOUBLE PRECISION, DIMENSION(nvar) :: phi

      INTEGER :: l,i


      Do i = 1,compteur_feuille     
! ============ Integration Radau5 ========================
       ITOLREAC = 0
       IJOMREAC = 0
       IMASREAC = 0
       MLMASREAC = 0
       MUMASREAC = 0
       IPARREAC = 0
       IOUTREAC = 0
       MLJACREAC = nvar
       MUJACREAC = nvar

       atol = 1.d-10
       rtol = 1.d-10
       !!atol = 1.d-15
       !!rtol = 1.d-15
 
       x = temps
       xend = x + dt/2.d0
       HRREAC = dt/2.d0
       RPARREAC = 0. 

      do l = 1,nvar
         phi(l) = tab_phi(i)%var(l)%t_ptr
      enddo
 

     CALL RADAU5(nvar,reaceq,x,phi,xend,HRREAC,rtol,&
                 atol,ITOLREAC,reaceq,IJOMREAC,MLJACREAC,&
                 MUJACREAC,reaceq,IMASREAC,MLMASREAC,MUMASREAC,&
                 reaceq,IOUTREAC,ELWRKREAC,LENWKREAC,IELWRKREAC,&
                 LENIWKREAC,RPARREAC,IPARREAC,ISTATEREAC)

      

      do l = 1,nvar
        if (phi(l) < 0d0) phi(l) = 0d0
        if (dabs(phi(l)-1d0) .le. zero) phi(l) = 1d0
        tab_phi(i)%var(l)%t_ptr = phi(l)
      enddo

      if (ISTATEREAC < 0) Print*,' u = ',phi 

      Enddo 

      end subroutine integration_reaction
