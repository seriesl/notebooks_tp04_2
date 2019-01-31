      subroutine lect_don

      use mod_common

      IMPLICIT NONE

      CHARACTER(len=80) :: nom_fich

      INTEGER :: l, lect

! +++++lecture des donnees dans le fichier don_euler.dat, unite=lect

      lect=1
      open(lect,file='don_burgers.dat',form='formatted')

! .....limites du domaine, nombre de racines
!      et flags des cond. aux limites

      nbre_racine(:) = 1

      do l = 1, ndim
        !!read(lect,*) lim_deb(l), lim_fin(l), nbre_racine(l)
        !!read(lect,*) lim_deb(l), lim_fin(l)
        lim_deb(l) = -2.
        lim_fin(l) = 2.
        !!read(lect,*) caltype_deb(l), caltype_fin(l)
        caltype_deb(l) = 2
        caltype_fin(l) = 2
      enddo

      !!Print*
      !!Print*, ' nombre de racines = ', nbre_racine

! .....initialisation de la solution a gauche
      !do l = 1, nvar
      !  read(lect,*) v_left(l)
      !enddo
      v_left(:) = 1.d0

! .....initialisation de la solution a droite
      do l = 1, nvar
      !!  read(lect,*) v_right(l)
        v_right(l) = 0.d0
      enddo
      

! .....Longueur de reference
!!      read(lect,*) L_ref
      L_ref = 1.

      if( L_ref == 0. ) L_ref = 1.

! .....Cle condition sur temperature de paroi
!!      read(lect,*) cle_twall
      cle_twall = .false.

! .....Temperature de paroi
!!      read(lect,*) t_wall
      t_wall = 1.
 
! .....Nombre de Mach
!!      read(lect,*) Mach
      Mach = 1.

! .....Nombre de Reynolds
!!      read(lect,*) Reynolds
      Reynolds = 0.

! .....Nombre de Prandtl
!!      read(lect,*) Prandtl
       Prandtl = 0.72

! .....nombre CFL
      !!read(lect,*) CFL
      CFL = 0.8


! .....choix du niveaux minimum de grilles
      read(lect,*) niv_min

! .....choix du nombre de niveaux de grilles
      read(lect,*) niv_max

! .....choix de l'ordre de reconstruction s 
      read(lect,*) s

! .....choix de la valeur d'epsilon 
      read(lect,*) epsilon

! .....nombre d'iteration en temps
      read(lect,*) npdt

! .....Cle de reprise ==> reprise = True
!!      read(lect,*) cle_reprise
      cle_reprise = .false.

! .....temps final de simulation
      read(lect,*) tps_fin

! .....premier temps de sauvegarde
      read(lect,*) tps_deb

! .....increment en temps pour la sauvegarde
      read(lect,*) dtinc
      !!print *, dtinc


      close(lect)

! +++++Valeurs de reference

      if( Reynolds == 0. ) then
        !Print*, ' Calcul Euler : ', ndim, ' D'
        !Print*, ' ------------ '
        !Print*
        !Print*, ' - Nombre de Mach      = ', Mach
      else
        Print*
        Print*
        Print*, ' Calcul Navier-Stokes : ', ndim, ' D'
        Print*, ' -------------------- '
        Print*
        Print*, ' - Nombre de Mach      = ', Mach
        Print*, ' - Nombre de Reynolds  = ', Reynolds
        Print*, ' - Nombre de Prandtl   = ', Prandtl
      endif

! .....nombre total de mailles

      nbre_maille = 2**(ndim * niv_max)
      Nbre_pt_arbre = ( 2**((niv_max+1)*ndim) - 1 ) / ( 2**ndim - 1 )
      
      do l = 1, ndim
        nbre_maille = nbre_maille * nbre_racine(l)
        Nbre_pt_arbre = Nbre_pt_arbre * nbre_racine(l)
      enddo

      print*
      print *, ' Multiresolution parameter:'
      print *, '   Minimal level:', niv_min
      print *, '   Maximal level:', niv_max
      print *, '   Predictor stencil :', s
      print *, '   Threshold value :', epsilon
      print *, '   Number of cells on the finest grid  : ', nbre_maille

      end subroutine lect_don
