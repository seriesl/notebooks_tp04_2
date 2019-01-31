      recursive subroutine construction_phi(maille,ll)

      use mod_common
      use mod_structure
      use mod_recherche

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      INTEGER, INTENT(INOUT) :: ll

      LOGICAL :: leave

      INTEGER :: i, l, m

! *****On marque les feuilles reelles (non-fictives)

      leave = .false.

      if( .not.associated( maille%fils(1)%ptr) ) then
        if( maille%feuille ) leave = .true.
      else
        if( .not.maille%fils(1)%ptr%feuille ) leave = .true.
      endif

! +++++On parcours les feuilles reelles (non-fictives)

      If( leave ) then

       ll = ll + 1
!......variables 

        do m = 1, ndim

! .....initialisation des flux

          do i = 1,nvar

             tab_phi(ll)%fdiff_dt(i,m)%t_ptr  =>  maille%fdiff_dt(i,m)
             tab_phi(ll)%fdiff_gch(i,m)%t_ptr =>  maille%fdiff_gch(i,m)

          enddo

!.......definition des voisins par dimension


          if( associated(maille%support(m)%precedent) ) then

                tab_phi(ll)%support(m)%precedent => maille%support(m)%precedent

          else
                  
                nullify(tab_phi(ll)%support(m)%precedent)

          endif

          if( associated(maille%support(m)%suivant) ) then

                tab_phi(ll)%support(m)%suivant => maille%support(m)%suivant

          else
                  
                nullify(tab_phi(ll)%support(m)%suivant)

          endif

          tab_phi(ll)%dx(m) = maille%dx(m)

        enddo

!...... definition des fils (feuilles fictives)

        if( associated(maille%fils(1)%ptr) ) then
 
          do m = 1, 2**ndim

            tab_phi(ll)%fils(m)%ptr => maille%fils(m)%ptr

          enddo

        else
       
          do m = 1, 2**ndim

            nullify(tab_phi(ll)%fils(m)%ptr)

          enddo


        endif


!...... definition des variables dans les feuilles

        
         do i = 1, nvar

            tab_phi( ll )%var(i)%t_ptr => maille%u(i)
     
         enddo
         
! +++++On remonte vers les fils

      Else

        if( associated(maille%fils(1)%ptr) ) then
          do l = 1, 2**ndim
            call construction_phi(maille%fils(l)%ptr, ll)
          enddo
        endif

      Endif

      end subroutine construction_phi
