      subroutine liberation(maille)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER :: l, m
      INTEGER :: indice, increment

      TYPE(structure_maille), pointer :: maille

      TYPE(structure_maille), pointer :: courant

! +++++On desalloue les fils et toute la structure associee

      do l = 1, 2**ndim

        if( associated(maille%fils(l)%ptr) ) then

          courant => maille%fils(l)%ptr

          deallocate(maille%fils(l)%ptr)

        endif

      enddo

      end subroutine liberation
