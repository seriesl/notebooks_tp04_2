      subroutine recherche_support(maille, s_st, f_support)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      TYPE(structure_maille), POINTER :: maille

      LOGICAL :: flag, info

      INTEGER :: ir, jr, kr
      INTEGER :: i, j, k, l, lm, m, n
      INTEGER :: i_sym, j_sym, k_sym
      INTEGER :: indice, increment
      INTEGER :: niv_sup, s_st, sj, sk
      INTEGER :: si_deb, sj_deb, sk_deb, si_fin, sj_fin, sk_fin
      INTEGER, DIMENSION(ndim) :: ind_loc, ind_sup, ind_cor
      INTEGER, DIMENSION(3) :: ind

      REAL(kind=8) :: rho_ec, rho_wall, rhoE_wall

      REAL(kind=8), DIMENSION(-s_st:s_st,-s_st:s_st,-s_st:s_st,nvar) :: &
                    f_support

      TYPE(structure_maille), pointer :: courant
      TYPE(structure_maille), pointer :: root_loc
 
! =======================================================================

! .....niveau et indice de la maille courante

      niv_sup = maille%maille_niv
      ind_sup(:) = maille%maille_i(:)

! +++++Definition des bornes debut et fin

        si_deb = max(1, ind_sup(1)-s_st) - ind_sup(1)
        si_fin = min(2**niv_sup*nbre_racine(1), ind_sup(1)+s_st) - ind_sup(1)

      if(ndim == 3) then
        sj = s_st
        sk = s_st

        sj_deb = max(1, ind_sup(2)-sj) - ind_sup(2)
        sj_fin = min(2**niv_sup*nbre_racine(2), ind_sup(2)+sj) - ind_sup(2)

        sk_deb = max(1, ind_sup(3)-sk) - ind_sup(3)
        sk_fin = min(2**niv_sup*nbre_racine(3), ind_sup(3)+sk) - ind_sup(3)

      elseif(ndim == 2) then
        sj = s_st
        sk = 0

        sj_deb = max(1, ind_sup(2)-sj) - ind_sup(2)
        sj_fin = min(2**niv_sup*nbre_racine(2), ind_sup(2)+sj) - ind_sup(2)

        sk_deb = 0
        sk_fin = 0

      elseif(ndim == 1) then
        sj = 0
        sk = 0

        sj_deb = 0
        sj_fin = 0

        sk_deb = 0
        sk_fin = 0

      endif

! ======== Recherche d'un support d'interpolation [-s_st , +s_st] ========

! +++++initialisation
      f_support(:,:,:,:) = 0.d0

! .....Valorisation pour la maille courante
      do m = 1,nvar
        f_support(0,0,0,m) = maille%u(m)
      enddo
          
! *****Recherche de tous les points du support et valorisation 

! .....polynome de degre >= 1
      if(s_st > 0 ) then

! +++++Recherche des points dan sle domaine de calcul

        do k = sk_deb, sk_fin
          ind(3) =  k

          do j = sj_deb, sj_fin
            ind(2) = j

            do i = si_deb, si_fin
              ind(1) = i

              do l = 1,ndim
                ind_loc(l) = ind_sup(l) + ind(l)
              enddo

! .....recherche de la racine appropriee
                  do kr = 1, nbre_racine(3)
                    do jr = 1, nbre_racine(2)
                      do ir = 1, nbre_racine(1)

                        flag = .true.

                        do l = 1, ndim
                         if( ind_loc(l) > &
              racine(ir,jr,kr)%ptr%maille_i(l)*2**niv_sup )  flag = .false.
                        enddo
 
                        if( flag ) goto 20

                      enddo
                    enddo
                  enddo

                  if( .not.flag ) then
                    Print*, ' Recherche Support : Probleme de ',&
                             'localisation racine !! '
                    Stop
                  endif

20                continue

! .....recherche du point dans la racine appropriee
                  root_loc => racine(ir,jr,kr)%ptr
                  call recherche(root_loc, ind_loc, niv_sup, courant, info)

              if( .not.info ) then

                print*, ' !! Recherche Support !! : ', &
                        ' Pb de recherche support maille - ', &
                        'Maille niv = ',niv_sup, ' Reelle = ',maille%feuille, &
                        ' Maille ijk = ',ind_sup, &
                        ' i, j, k loc = ',ind_loc

                stop
              endif


              do m = 1,nvar
                f_support(i, j, k, m) = courant%u(m)
              enddo

            enddo
          enddo
        enddo

! ========================= CAL suivant X ==========================

! +++++Conditions aux limites amont

        if( si_deb > -s_st ) then


            do k = -sk, sk
              do j = -sj, sj

                do i = -s_st, si_deb-1
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(si_deb,j,k,m)
                  enddo
                enddo

              enddo
            enddo


        endif

! *****Conditions a la limite aval

        if( si_fin < s_st ) then

            do k = -sk, sk
              do j = -sj, sj

                do i = si_fin+1, s_st
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(si_fin,j,k,m)
                  enddo
                enddo

              enddo
            enddo


        endif

! ========================= CAL suivant Y ==========================

! +++++Conditions aux limites amont

        if( sj_deb > -sj ) then

            do k = -sk, sk
              do i = -s_st, s_st

                do j = -sj, sj_deb-1
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,sj_deb,k,m)
                  enddo
                enddo

              enddo
            enddo

        endif

! *****Conditions a la limite aval

        if( sj_fin < sj ) then

            do k = -sk, sk
              do i = -s_st, s_st

                do j = sj_fin+1, sj
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,sj_fin,k,m)
                  enddo
                enddo

              enddo
            enddo

        endif

! ========================= CAL suivant Z ==========================

! +++++Conditions aux limites amont
        if( sk_deb > -sk ) then

            do j = -sj, sj
              do i = -s_st, s_st

                do k = -sk, sk_deb-1
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,j,sk_deb,m)
                  enddo
                enddo

              enddo
            enddo

       endif
! *****Conditions a la limite aval

        if( sk_fin < sk ) then

! .....Conditions de frontiere fluide ==> derivee nulle
!          elseif( caltype_fin(3) == 0 ) then

            do j = -sj, sj
              do i = -s_st, s_st

                do k = sk_fin+1, sk
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,j,sk_fin,m)
                  enddo
                enddo

              enddo
            enddo

       endif

      endif
      end subroutine recherche_support
