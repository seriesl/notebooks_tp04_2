      recursive subroutine seuillage(maille)

      use mod_common
      use mod_structure
      use mod_recherche

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      LOGICAL :: flag, info

      INTEGER :: ir, jr, kr
      INTEGER :: i, j, k, l, m
      INTEGER :: indice, increment
      INTEGER :: niv_sup, sj, sk
      INTEGER, DIMENSION(ndim) :: ind_loc, ind_sup, ind_cor
      INTEGER, DIMENSION(3) :: ind

      TYPE(structure_maille), pointer :: courant
      TYPE(structure_maille), pointer :: root_loc
 
      REAL(KIND = 8) :: precis, eps

! ============= seuillage et marquage des branches a elaguer =========

! *****C'est une feuille

      if( .not.associated(maille%fils(1)%ptr) ) then

! .....des fils ne doivent pas etre crees

        maille%arbre = .false.

! .....Mais si le pere d'une feuille reelle a un detail important 
!           ==> creation potentielle

        if( maille%feuille ) then        

! +++++Valorisation de seuil suivant le niveau du pere
! L1 norm
!!        eps = epsilon * ( 2.**(ndim*((maille%maille_niv-1) - niv_max)) )
 
! L2 norm
!        eps = epsilon * dsqrt( 2d0**dble(ndim*((maille%maille_niv-1) - niv_max))) ! detail => pere
       eps = epsilon * dsqrt( 2d0**dble(ndim*(maille%maille_niv - niv_max))) ! detail => feuille moy

 

!          if((abs(maille%pere%detail / det_max) >= eps) ) &
          if(dabs(maille%pere%detail)>=eps ) &
             maille%arbre = .true.

        endif
        
      else

! *****On passe aux fils

        do l = 1, 2**ndim
         call seuillage(maille%fils(l)%ptr)
        enddo

! *****test sur le pere

! .....le detail est negligeable ==> point potentiellement plus dans arbre

! +++++Valorisation de seuil suivant le niveau courant
! L1 norm
!!        eps = epsilon * ( 2.**(ndim*(maille%maille_niv - niv_max)) )

! L2 norm
!        eps = epsilon * dsqrt( 2d0**dble(ndim*(maille%maille_niv - niv_max))) ! detail => pere
        eps = epsilon * dsqrt( 2d0**dble(ndim*((maille%maille_niv+1) - niv_max))) ! detail => feuille moy


!        if((abs(maille%detail / det_max) < eps) &
        if(dabs(maille%detail) < eps  &
         .and. (maille%maille_niv > niv_min)) then

          maille%arbre = .false.
!!!!!!!!!!!!!detail => feuille moy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do l = 1, 2**ndim
             maille%fils(l)%ptr%arbre = .false.
            enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ONLY IF THE LEAF IS REAL
!          if (.not.maille%fils(1)%ptr%feuille) print*,'det =',maille%detail 

! .....sinon, le point est potentiellement dans arbre
        else
          maille%arbre = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! ALTERNATIVE TO HARTEN APPORACH => SECURITY ZORE
!!!!!!!!!! WE ADD ANOTHER LAYER
!!!!!!!!!!!!!detail => feuille moy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do l = 1, 2**ndim
             maille%fils(l)%ptr%arbre = .true.
            enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! .....niveau et indice de la maille courante
          niv_sup = maille%maille_niv
          ind_sup(:) = maille%maille_i(:)


! .....ainsi que ces proches voisins pour reconstruction polynomiale
!      si polynome de degre >= 1
          if(s > 0 ) then

            if(ndim == 3) then
              sj = s
              sk = s
            elseif(ndim == 2) then
              sj = s
              sk = 0
            elseif(ndim == 1) then
              sj = 0
              sk = 0
            endif

! .....Recherche de tous les points du support
            do k = -sk, sk
              ind(3) =  k
              do j = -sj, sj
                ind(2) = j
                do i = -s, s
                  ind(1) = i

                  do l = 1,ndim
                    ind_loc(l) = ind_sup(l) + ind(l)
                  enddo

! .....conditions aux limites
                  call cal_indice(niv_sup, ind_loc, ind_cor)

! .....recherche de la racine appropriee
                  do kr = 1, nbre_racine(3)
                    do jr = 1, nbre_racine(2)
                      do ir = 1, nbre_racine(1)

                        flag = .true.

                        do l = 1, ndim
                         if( ind_cor(l) > &
              racine(ir,jr,kr)%ptr%maille_i(l)*2**niv_sup )  flag = .false.
                        enddo
 
                        if( flag ) goto 20

                      enddo
                    enddo
                  enddo

                  if( .not.flag ) then
                    Print*, ' Seuillage : Probleme de localisation ',&
                             'racine !! '
                    Stop
                  endif

20                continue

! .....recherche du point dans la racine appropriee
                  root_loc => racine(ir,jr,kr)%ptr
                  call recherche(root_loc, ind_cor, niv_sup, courant, info)

! +++++le point est trouve dans l'arbre
                  if( info ) then

                    courant%arbre = .true.

! +++++le point n'est pas dans l'arbre !!
                  else
                    print*, ' !! Seuillage !! : ', &
                            ' Pb de recherche support maille - ', &
                            'Maille niv = ',niv_sup,' Maille ijk = ',ind_sup, &
                            ' i, j, k loc = ',ind_loc,' i, j, k cor = ',ind_cor
                    stop
                  endif



                enddo
              enddo
            enddo

          endif

! +++++test sur croissance des singularites

! +++++Valorisation de seuil suivant le niveau courant
  
!!        eps = epsilon * ( 2.**(ndim*(maille%maille_niv - niv_max)) )


!!!!!!!!!!!!!!!!!HARTEN APPORACH WITH NIV+2!!!!!!!!!!!!!!!!!!!!!
! L2 norm
!       eps = epsilon * dsqrt( 2d0**dble(ndim*(maille%maille_niv - niv_max))) ! detail => pere
!        eps = epsilon * dsqrt( 2d0**dble(ndim*((maille%maille_niv+2) - niv_max))) ! detail => feuille moy


!----------norm L1
!!          precis = 2.**(2*s-1) * eps

!----------norm L2
!          detail => pere
!!!!!!!!!!!!HARTEN APPROACH
!           precis = 2.d0**(2.d0*dble(s)) * eps !+dble(ndim)/2.d0) * eps !harten p+1=2s+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          precis = 2.d0**(dble(2*s-1)/2.d0) * eps !CT

!          detail => feuille moy ADD PRECIS FOR MORE COMPRESSION INSTEAD
!          OF BOUCLE BEFORE
!           precis = dsqrt(2.d0**dble(ndim)) * eps ! threshold of next level harten p=-1


!!          if(abs(maille%detail / det_max) >= precis ) then
!          if(dabs(maille%detail) >= precis ) then

!            maille%arbre = .true.

!            do l = 1, 2**ndim
!              maille%fils(l)%ptr%arbre = .true.
!            enddo

!          endif

        endif  

      endif

      end subroutine seuillage
