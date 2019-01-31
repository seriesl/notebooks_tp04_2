      subroutine creation(maille, reelle)

      use mod_common
      use mod_structure

      IMPLICIT NONE	

      LOGICAL :: reelle

      TYPE(structure_maille), pointer :: maille

      INTEGER :: i, k, l, lm, m, n, ll, mm, nn, kk
      INTEGER :: indice, increment

      DOUBLE PRECISION :: delta

! +++++On alloue et on valorise les fils

      do l = 1, 2**ndim

! *****Allocation de maille_fils

        allocate(maille%fils(l)%ptr)

! *****valorisation des valeurs de la structure

! .....niveau
        maille%fils(l)%ptr%maille_niv= maille%maille_niv+1		

! .....indices
        k = 0
        n = 0
        do m = ndim, 1, -1
          k = k+1
          i = ((l-1) - n)/2**(ndim-k)
          n = n + i*2**(ndim-k)

          maille%fils(l)%ptr%maille_i(m) = 2*maille%maille_i(m)+i-1

! .....coordonnees
          delta = ( dble(i)*.5d0 - 0.25d0 ) * maille%dx(m)
          maille%fils(l)%ptr%x(m)= maille%x(m) + delta

! .....pas d'espace
          maille%fils(l)%ptr%dx(m)= .5d0*maille%dx(m)

        enddo
        
! .....la maille a-t-elle des fils ?
        maille%fils(l)%ptr%arbre = .false.
        
! .....la feuille est-elle reelle ou fictive
        maille%fils(l)%ptr%feuille = reelle

! .....details
        do lm = 1, 2**ndim
          do m = 1,nvar
            maille%fils(l)%ptr%det(m, lm) = 0.d0
          enddo
        enddo

        maille%fils(l)%ptr%detail = 0.d0

! .....initialisation des flux

        maille%fils(l)%ptr%fdiff_dt(:,:) = 0.d0
        maille%fils(l)%ptr%fdiff_gch(:,:) = 0.d0
                
! *****chainage	

        maille%fils(l)%ptr%pere => maille		

      enddo

      end subroutine creation
