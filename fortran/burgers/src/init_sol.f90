      recursive subroutine init_sol(maille)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      EXTERNAL cellave2

      INTEGER :: l, m
      REAL(kind=8) :: x
      REAL(kind=8) :: rho_ec, gh_inf
      REAL(kind=8) :: rhot_inf, pt_inf, p_loc, s_choc
      REAL(kind=8) :: R_0, rayon, ampli, c2, u_theta, t_loc, delta_temp
      REAL(kind=8) :: rho_aval, p_aval
      REAL(kind=8), DIMENSION(ndim) :: x_0, x_low, x_high

      TYPE(structure_maille), pointer:: maille

! *****traitement des feuilles

      if(.not.associated(maille%fils(1)%ptr)) then

! +++++Initialisation globale

        if( ndim == 1 ) then
          do m = 1, ndim
            !!x_0(m) = 0.5*( lim_fin(m) + lim_deb(m) )
            x = maille%x(m)
            do l = 1, nvar
              !maille%u(l) = - v_left(l)*sin( 2.*pi*(maille%x(m)-x_0(m))/ &
              !                               (lim_fin(m) - lim_deb(m)) ) 
              if (x > -1.d0 .and. x <= 0.d0) then
                maille%u(l) = maille%x(m) + 1.d0
              else if (x > 0.d0 .and. x < 1.d0) then 
                maille%u(l) = 1.d0 - maille%x(m) 
              else  
                maille%u(l) = 0.d0
              endif  
            enddo
            !!print *, maille%x(m), maille%u(1)
          enddo

        else

! .....fraction de l'etat gauche c2
          do m = 1, ndim
            x_low(m) = maille%x(m) - maille%dx(m)*.5
          enddo
          call cellave2(x_low(1), x_low(2), maille%dx(1), maille%dx(2),c2)

          t_loc = (1.-c2)*1. + c2*1000./300.
          do l = 1, nvar
            maille%u(l) = (1.-c2)*v_left(l) + c2*v_right(l)
          enddo

        endif

      else

! +++++traitement des fils

        if( associated(maille%fils(1)%ptr) ) then
          do l = 1, 2**ndim
            !!print *, '  l =', l 
            call init_sol(maille%fils(l)%ptr)
          enddo
        endif

      endif

      end subroutine init_sol

