      recursive subroutine built_fine_sol(maille)

      use mod_common
      use mod_structure
      use mod_recherche

      IMPLICIT NONE
   
      LOGICAL :: flag, reelle

      INTEGER :: l, m

      TYPE(structure_maille), pointer :: maille

! +++++On marque les feuilles reelles (non-fictives)

      flag = .false.

      if( .not.associated( maille%fils(1)%ptr) ) then
        if( maille%feuille ) flag = .true.
      else
        if( .not.maille%fils(1)%ptr%feuille ) flag = .true.
      endif

! +++++On cree les enfants des feuilles reelles (non-fictives)

      If( flag .and. maille%maille_niv < niv_max) then

          reelle = .true.
         if( .not.associated( maille%fils(1)%ptr) )  then
              call creation(maille, reelle)
         else
              do l = 1, 2**ndim
                   maille%fils(l)%ptr%feuille = reelle
              enddo
         endif

         do l = 1, 2**ndim
            do m = 1,nvar
               maille%det(m, l) = 0d0
            enddo
         enddo

          compteur_create = compteur_create + 1

          call graduation_local(maille, reelle)

! .....Valorisation des Fils
          call decoding_moy(maille)


      Else

! +++++On remonte vers les fils

        if( associated(maille%fils(1)%ptr) ) then
          do l = 1, 2**ndim
            call built_fine_sol(maille%fils(l)%ptr)
          enddo
        endif

      Endif


      end subroutine built_fine_sol
