      recursive subroutine compute_error(maille)

      use mod_common
      use mod_structure

      IMPLICIT NONE
   
      INTEGER :: l, m

      TYPE(structure_maille), pointer :: maille

! +++++Calcule erreur sur la grille la plus fine

      If( .not.associated( maille%fils(1)%ptr) ) then

        if (save_res == 1 .or. dabs(temps-tps_fin).le.zero) then 
         write(32,102) ( maille%x(m),m=1,ndim ), &
           ( maille%u(m)-sol_split(nvar*i_sav+m,i_err),m=1,nvar ), &
           ( maille%u(m),m=1,nvar),(sol_split(nvar*i_sav+m,i_err),m=1,nvar )
        endif

        do m = 1,nvar
           err_max(m) = &
            max( err_max(m),dabs(sol_split(nvar*i_sav+m,i_err)-maille%u(m)))! &
!             / sol_split(nvar*i_sav+m,i_err)) 

           error_l2(m) = error_l2(m) + &
               (sol_split(nvar*i_sav+m,i_err)-maille%u(m)) * &
               (sol_split(nvar*i_sav+m,i_err)-maille%u(m)) !&
!            / (sol_split(nvar*i_sav+m,i_err) * &
!               (sol_split(nvar*i_sav+m,i_err))))
         enddo

         i_sav = i_sav + 1

      Else

! +++++On remonte vers les fils

        if( associated(maille%fils(1)%ptr) ) then
          do l = 1, 2**ndim
            call compute_error(maille%fils(l)%ptr)
          enddo
        endif

      Endif

!+++for 3D problems    
!102   format(3(1x,1pe15.8),1x,i4,6(1x,1pe15.8))

!+++for 2D problems    
!102   format(2(1x,1pe15.8),1x,i4,5(1x,1pe15.8))

!+++for 1D problems    
102   format(4(1x,1pe15.8))

      end subroutine compute_error
