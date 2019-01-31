      subroutine print_leaves

      use mod_common
      use mod_structure
      use mod_fonction

      IMPLICIT NONE
   
      INTEGER :: i, ijk, j, k, l, leaf, m, n
      INTEGER :: Jmj
      INTEGER, DIMENSION(3) :: i_c
      INTEGER, DIMENSION(8,3) :: ordre_coins

      REAL(kind=8) :: rayon, ampli, u_theta, t_loc, delta_temp
      REAL(kind=8) :: rho_exact
      REAL(kind=8), DIMENSION(ndim) :: x_0, u_exact

      TYPE(structure_maille), pointer :: maille
      REAL(kind=8) :: bord_g, bord_d

! +++++definition de l'ordre de parcours des coins
      ordre_coins(1,1) = 0
      ordre_coins(1,2) = 0
      ordre_coins(1,3) = 0
      ordre_coins(2,1) = 1
      ordre_coins(2,2) = 0
      ordre_coins(2,3) = 0
      ordre_coins(3,1) = 1
      ordre_coins(3,2) = 1
      ordre_coins(3,3) = 0
      ordre_coins(4,1) = 0
      ordre_coins(4,2) = 1
      ordre_coins(4,3) = 0 
      ordre_coins(5,1) = 0
      ordre_coins(5,2) = 0
      ordre_coins(5,3) = 1
      ordre_coins(6,1) = 1
      ordre_coins(6,2) = 0
      ordre_coins(6,3) = 1
      ordre_coins(7,1) = 1
      ordre_coins(7,2) = 1
      ordre_coins(7,3) = 1
      ordre_coins(8,1) = 0
      ordre_coins(8,2) = 1
      ordre_coins(8,3) = 1 

! +++++Eciture au format colonne pour le 1D

      if( ndim == 1 ) then

! *****On parcourt les feuilles reelles (non-fictives)
        Do leaf = 1, compteur_feuille

          maille => feuille_reelle(leaf)%ptr
!!
          if( maille%feuille ) then
!!
! +++++Ecriture des coordonnees de l'arbre et des variables
            do m = 1, ndim
              x_coins(leaf, m) = maille%x(m)
            enddo

            var_print(leaf, 1) = maille%dx(1)
            var_print(leaf, 2) = maille%maille_niv
            do m = 1, nvar
              var_print(leaf, m+2) = maille%u(m)
              !!print *, var_print(leaf, m+2)  
            enddo

            !!if (temps == 0.d0) then
            !!  var_print(leaf, 4) = 0.d0
            !!else if (temps <= 1.d0) then
            !!  bord_g = maille%x(1)-(0.5*maille%dx(1)) 
            !!  bord_d = maille%x(1)+(0.5*maille%dx(1)) 
            !!  if (maille%x(1) < -1.d0  ) then
            !!    print *, "LS: err0: ", maille%x(1), 0.d0, maille%u(1)
            !!    var_print(leaf, 4) = dabs(maille%u(1) - 0.d0)
            !!  !!else if (bord_g > -1.d0 .and. bord_d <= temps) then
            !!  else if (maille%x(1) > -1.d0 .and. maille%x(1) < temps) then
            !!    print *, "LS: err1: ", maille%x(1), (maille%x(1)+1.d0)/(1.d0+temps), maille%u(1)
            !!    var_print(leaf, 4) = dabs(maille%u(1) - ((maille%x(1)+1.d0)/(1.d0+temps)))
            !!  !!else if ((bord_g >= temps .and. bord_d<=temps) .or. (bord_g==temps)) then
            !!  !!else if ((maille%x(1) >= temps .and. bord_d<=temps) .or. (bord_g==temps)) then
            !!    !!print *, "LS: err5: ", maille%x(1), (maille%x(1)+1.d0)/(1.d0+temps), maille%u(1)
            !!    !!var_print(leaf, 4) = dabs(maille%u(1) - 1.d0)
            !!  else if (maille%x(1) > temps .and. maille%x(1) < 1.d0) then
            !!    print *, "LS: err2: ", maille%x(1),(1.d0-maille%x(1))/(1.d0-temps), maille%u(1)
            !!    var_print(leaf, 4) = dabs(maille%u(1) - ((-maille%x(1)+1.d0)/(1.d0-temps)))
            !!  else   
            !!    print *, "LS: err3: ",  maille%x(1), 0.d0, maille%u(1)
            !!    var_print(leaf, 4) = dabs(maille%u(1) - 0.d0)
            !!  end if

            !!endif
!!
          endif
!!
        Enddo
    
! +++++Ecriture au format E-F pour le Multi-D       

      Else

! *****On parcourt les feuilles reelles (non-fictives)

        Do leaf = 1, compteur_feuille

          maille => feuille_reelle(leaf)%ptr
!!
          if( maille%feuille ) then
!!
! *****Ecriture des coordonnees de l'arbre et des variables

            Jmj = niv_max - maille%maille_niv

! *****Pour les 2**ndim coins du volume de controle
            do l = 1, 2**ndim

! .....calcul du mono-indice
              i_c(:) = 0
              do m = 1, ndim
                i_c(m) = 2**Jmj*maille%maille_i(m) + &
                         (ordre_coins(l,m)-1)*2**Jmj
              enddo
              ijk = i_c(1)+1+i_c(2)*(i_max+1)+i_c(3)*(i_max+1)*(j_max+1)

! .....Ecriture des elements
              if( flag_coins(ijk) < 0 ) then
                nbre_coins = nbre_coins + 1
                flag_coins(ijk) = nbre_coins
                do m = 1, ndim
                  x_coins(nbre_coins, m) = maille%x(m) + &
                    (real(ordre_coins(l,m))-0.5)*maille%dx(m)
                enddo
              endif
              ind_coins(leaf, l) = flag_coins(ijk)

            enddo

! .....Ecriture des variables
            do m = 1, nvar
              var_print(leaf, m) = maille%u(m)
            enddo

!!
          endif
!!
        Enddo

      Endif

      end subroutine print_leaves
