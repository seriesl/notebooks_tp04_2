      recursive subroutine init_flags(maille)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER :: l, indice, increment
      TYPE(structure_maille), pointer :: maille


! +++++Initialisation de la table de hashage

      If ( associated(maille) ) then

! +++++Initialisation des mailles dans l'arbre

        maille%arbre = .false.

! .....On descend vers les fils
        if(associated(maille%fils(1)%ptr)) then

          do l = 1, 2**ndim
            call init_flags(maille%fils(l)%ptr)
         enddo

        endif

      Endif

      end subroutine init_flags
