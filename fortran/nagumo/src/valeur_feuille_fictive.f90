      subroutine valeur_feuille_fictive(mm)

      use mod_common
      use mod_structure
      use mod_recherche

      implicit none


      INTEGER :: ll, feu, mm ! var separees

      TYPE(structure_maille), pointer :: maille


      INTEGER :: i, j, k, l, m, n
      INTEGER :: indice, increment
      INTEGER :: niv_sup, sj, sk
      INTEGER, DIMENSION(ndim) :: ind_loc, ind_sup, ind_cor
      INTEGER, DIMENSION(3) :: ind

      DOUBLE PRECISION, DIMENSION(nvar) :: qx, qy, qz,  qxy, qyz, qxz, qxyz
      DOUBLE PRECISION, DIMENSION(nvar) :: u_tilde

      DOUBLE PRECISION,DIMENSION(-s:s,-s:s,-s:s,nvar) :: f_support


! +++++On valorise les feuilles fictives a partir des peres

!......On parcours tab_phi

       Do ll = 1,compteur_feuille
     
       If( associated( tab_phi(ll)%fils(1)%ptr) ) then


          maille =>  tab_phi(ll)%fils(1)%ptr%pere


!......valeurs des feuilles fictives obtenues par interpolation
! *****Valorisation des fils fictifs


! .....Allocation du support
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


! .....niveau et indice de la maille courante
        niv_sup = maille%maille_niv
        ind_sup(:) = maille%maille_i(:)

!.....recherche du support de l'interpolation 

        if( s > 0 ) call recherche_support_diff(maille, s, f_support, mm)



! *****calcul des polynomes d'interpolation

        if(s > 0) then
          do m = 1,nvar
            qx(m) = 0.d0

            do i = 1, s
              qx(m) = qx(m) + coef(i)* &
                      (f_support(i,0,0,m)-f_support(-i,0,0,m))
            enddo
          enddo
        else
          do m = 1,nvar
            qx(m) = 0.d0
          enddo
        endif

        if(sj > 0 ) then
          do m = 1,nvar
            qy(m) = 0.d0
            qxy(m) = 0.d0

            do j = 1, sj
              qy(m) = qy(m) + coef(j)* &
                      (f_support(0,j,0,m)-f_support(0,-j,0,m))

              do i = 1,s
                qxy(m) = qxy(m) + coef(j)*coef(i)* &
                         (f_support(i,j,0,m)-f_support(i,-j,0,m) - &
                          f_support(-i,j,0,m)+f_support(-i,-j,0,m))
              enddo
            enddo
          enddo
        else
          do m = 1,nvar
            qy(m) = 0.d0
            qxy(m) = 0.d0
          enddo
        endif

        if(sk > 0 ) then
          do m = 1,nvar
            qz(m) = 0.d0
            qxz(m) = 0.d0
            qyz(m) = 0.d0
            qxyz(m) = 0.d0

            do k = 1, sk
              qz(m) = qz(m) + coef(k)* &
                      (f_support(0,0,k,m)-f_support(0,0,k,m))

              do i = 1,s
                qxz(m) = qxz(m) + coef(k)*coef(i)* &
                         (f_support(i,0,k,m)-f_support(i,0,-k,m) - &
                          f_support(-i,0,k,m)+f_support(-i,0,-k,m))
              enddo

              do j = 1,sj
                qyz(m) = qyz(m) + coef(j)*coef(k)* &
                         (f_support(0,j,k,m)-f_support(0,-j,k,m) - &
                          f_support(0,j,-k,m)+f_support(0,-j,-k,m))

                do i = 1,s
                  qxyz(m) = qxyz(m) + coef(k)*coef(j)*coef(i)* &
                           (f_support(i,j,k,m)-f_support(i,j,-k,m) - &
                            f_support(i,-j,k,m)-f_support(-i,j,k,m) + &
                            f_support(i,-j,-k,m)+f_support(-i,j,-k,m) + &
                            f_support(-i,-j,k,m)-f_support(-i,-j,-k,m))
                enddo
              enddo
            enddo
          enddo
        else
          do m = 1,nvar
            qz(m) = 0.d0
            qxz(m) = 0.d0
            qyz(m) = 0.d0
            qxyz(m) = 0.d0
          enddo
        endif
            
! *****calcul de utilde par polynome d'interpolation

        do l = 1, 2**ndim
            
          ind(:) = 0
          k = 0
          n = 0
          do m = ndim, 1, -1
            k = k+1
            i = ((l-1) - n)/2**(ndim-k)
            n = n + i*2**(ndim-k)

            ind(m) = i
          enddo

          do m = 1,nvar
            u_tilde(m) = f_support(0, 0, 0, m)

! .....prediction sur le maillage grossier

            u_tilde(m) = u_tilde(m) - ( &
             - dble((-1)**ind(1))*qx(m) - dble((-1)**ind(2))*qy(m) - dble((-1)**ind(3))*qz(m) &
             + dble((-1)**(ind(1)+ind(2)))*qxy(m) + dble((-1)**(ind(2)+ind(3)))*qyz(m) &
             + dble((-1)**(ind(2)+ind(3)))*qxz(m) &
             - dble((-1)**(ind(1)+ind(2)+ind(3)))*qxyz(m) )

! .....valeurs sur les fils (au niveau plus fin)

            maille%fils(l)%ptr%u(m) = u_tilde(m)

          enddo

        enddo

  
      Endif

      Enddo

      end subroutine valeur_feuille_fictive
