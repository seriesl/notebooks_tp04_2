      recursive subroutine init_sol(maille)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER :: l, m, n

      DOUBLE PRECISION :: ko, beta, delta !KPP

      TYPE(structure_maille), pointer:: maille

! *****traitement des feuilles

      if(.not.associated(maille%fils(1)%ptr)) then

! +++++Initialisation sur la grille la plus fine

!!!!!! KPP 1D
       delta =  lim_deb(1) + (lim_fin(1)-lim_deb(1))/3.d0 

       if (kpp_ini == 0) then

         ko = 1./Dif(1)

         beta = min(1.d0,dexp(-(maille%x(1)-delta)*dsqrt(ko/(Dif(1)*2.d0)))/ &
           (1.d0+dexp(-(maille%x(1)-delta)*dsqrt(ko/(Dif(1)*2.d0)))))

         maille%u(1) = max(beta,zero)

         maille%u(1) = min(1d0,beta)

         if (dabs(maille%u(1) - 1d0) .le. zero ) maille%u(1) = 1d0
       else

         if (maille%x(1) .le. delta ) then 
             maille%u(1) = 1.d0 
         else
             maille%u(1) = 0.d0
         endif

       endif

 
      else

! +++++traitement des fils

        if( associated(maille%fils(1)%ptr) ) then
          do l = 1, 2**ndim
            call init_sol(maille%fils(l)%ptr)
          enddo
        endif

      endif

      end subroutine init_sol

