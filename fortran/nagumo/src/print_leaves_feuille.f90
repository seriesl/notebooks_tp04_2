      recursive subroutine print_leaves(maille)

      use mod_common
      use mod_structure

      IMPLICIT NONE
   
      LOGICAL :: flag

      INTEGER :: l, m

      DOUBLE PRECISION :: rho_ec, PsgM2, TsgM2
      DOUBLE PRECISION :: vorticite
      DOUBLE PRECISION, DIMENSION(ndim) :: grad_rho, grad_P, grad_T
      DOUBLE PRECISION, DIMENSION(ndim,ndim) :: grad_u

      DOUBLE PRECISION :: rayon, ampli, u_theta, t_loc, delta_temp
      DOUBLE PRECISION :: rho_exact
      DOUBLE PRECISION, DIMENSION(ndim) :: x_0, u_exact

      TYPE(structure_maille), pointer :: maille

! +++++On marque les feuilles reelles (non-fictives)

      flag = .false.

      if( .not.associated( maille%fils(1)%ptr) ) then
        if( maille%feuille ) flag = .true.
      else
        if( .not.maille%fils(1)%ptr%feuille ) flag = .true.
      endif

! +++++On ecrit les feuilles reelles (non-fictives)

      If( flag ) then

       if (save_res == 1 .or. dabs(temps-tps_fin).le.zero) then
          write(31,101) ( maille%x(m),m=1,ndim ), maille%maille_niv, &
                          maille%detail, & ! detail => feuille
!                          maille%pere%detail/2d0**dble(ndim), & ! detail => pere
!                        ( maille%dx(m),m=1,ndim ), & 
                        ( maille%u(m),m=1,nvar )
       endif


          if (i_err .gt. 0 .and. dabs(epsilon-0D0).le. zero) then
            do m=1,nvar
             sol_split(nvar*i_sav+m,i_err) = maille%u(m)
            enddo
            i_sav = i_sav + 1
          endif


      Else

! +++++On remonte vers les fils

        if( associated(maille%fils(1)%ptr) ) then
          do l = 1, 2**ndim
            call print_leaves(maille%fils(l)%ptr)
          enddo
        endif

      Endif

!+++for 3D problems    
!101   format(3(1x,1pe15.8),1x,i4,6(1x,1pe15.8))

!+++for 2D problems    
!101   format(2(1x,1pe15.8),1x,i4,5(1x,1pe15.8))

!+++for 1D problems    
101   format(1(1x,1pe15.8),1x,i4,4(1x,1pe15.8))

102   format(1x,2(1x,1pe15.8),1x,i4,7(1x,1pe15.8))

103   format(1x,10(1x,1pe15.8))

105   format(1x,5(1x,1pe15.8))

      end subroutine print_leaves
