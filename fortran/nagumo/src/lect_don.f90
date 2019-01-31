      subroutine lect_don

      use mod_common

      IMPLICIT NONE

      CHARACTER(len=80) :: nom_fich

      INTEGER :: l, lect

! +++++lecture des donnees dans le fichier don_euler.dat, unite=lect

      lect=1
      open(lect,file='don_nagumo.dat',form='formatted')


! .....limites du domaine, nombre de racines
!      et flags des cond. aux limites

      nbre_racine(:) = 1

      do l = 1, ndim
        !read(lect,*) lim_deb(l), lim_fin(l), nbre_racine(l)
        lim_deb(l) = -70.
        lim_fin(l) = 70.
        
!        lim_deb(l) = lim_deb(l) - (size/Lll)*0.5
!        lim_fin(l) = lim_fin(l)*(size/Lll)*0.5

      enddo

      !Print*
      !Print*, ' number of roots : ', nbre_racine
      
      !Print*
      !Print*, ' left boundary : ', lim_deb, ' right boundary = ', lim_fin



! .....choix du niveaux minimum de grilles
      read(lect,*) niv_min

! .....choix du nombre de niveaux de grilles
      read(lect,*) niv_max

! .....choix de l'ordre de reconstruction s 
      read(lect,*) s

! .....choix de la valeur d'epsilon 
!      read(lect,*) epsilon

      read(lect,*) eps_max
      !!read(lect,*) eps_min
      eps_min = eps_max

! .....nombre d'iteration en temps
      read(lect,*) npdt

! .....temps final de simulation
      read(lect,*) tps_fin

! .....premier temps de sauvegarde
      read(lect,*) tps_deb

! .....increment en temps pour la sauvegarde
      read(lect,*) dtinc

!.....save_res : save 0:only last result - 1: all results for t_delta
      !!read(lect,*) save_res
      save_res = 1


!......coefficients de diffusion (a, b, c)/(Ox, Fu, T)
      do l = 1, nvar
        read(lect,*) Dif(l)
      enddo

!=======KPP========================================
!.....initial condition

      !!read(lect,*) kpp_ini
      kpp_ini = 0

      close(lect)




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
      print *, '   Threshold value :', eps_max
      print *, '   Number of cells on the finest grid  : ', nbre_maille
      print *

      end subroutine lect_don
