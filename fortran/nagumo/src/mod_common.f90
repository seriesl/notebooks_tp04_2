!	module de declaration des variables communes

      module mod_common

      implicit none

! .....dimension spatiale du probleme
      INTEGER, PARAMETER :: ndim = 1
! .....nombre de variables
      INTEGER, PARAMETER :: nvar = 1

! .....ordre du scheme p = 2*nordre-1
!      INTEGER, PARAMETER :: nordre = 4
! .....niveau maximum des grilles emboitees
      INTEGER, PARAMETER :: niveau_maximum = 12

! .....parametre de precision
      DOUBLE PRECISION, PARAMETER :: zero = 1.d-14

!====== Parametres BZ=================================

!      DOUBLE PRECISION, PARAMETER :: k1 = 1.e+2 ! 1 / epsilon
!      DOUBLE PRECISION, PARAMETER :: k2 = 1.e+5 ! 1 / mhu
!      DOUBLE PRECISION, PARAMETER :: q = 2.e-3 ! q
!      DOUBLE PRECISION, PARAMETER :: ff = 16.e-1 ! f

!      DOUBLE PRECISION, DIMENSION(nvar) :: Dif  
!------------------------------ coefficient de diffusion (a, b, c) 

!====== Parametres Chimie Simple==============================

!      DOUBLE PRECISION, PARAMETER :: Bb = 1.e+7
!      DOUBLE PRECISION, PARAMETER :: Ta = 10055.
!      DOUBLE PRECISION, PARAMETER :: phii = 8.e-1
!      DOUBLE PRECISION, PARAMETER :: Tinit = 300.
!      DOUBLE PRECISION, PARAMETER :: Qq = 34550.
!      DOUBLE PRECISION, PARAMETER :: Kappa = 2.26d-5
!      DOUBLE PRECISION, PARAMETER :: Lll = 4.e-3
!      DOUBLE PRECISION, PARAMETER :: Lf = 1.
!      DOUBLE PRECISION, PARAMETER :: pm = 35d-2
!      DOUBLE PRECISION, PARAMETER :: size = 5d-2
!      DOUBLE PRECISION :: cst1, cst2, sizex, sizey 



!      DOUBLE PRECISION, DIMENSION(nvar) :: Dif  
!------------------------------ coefficient de diffusion (Ox, Fu, T) 

!====== Parametres KPP================================

      DOUBLE PRECISION, DIMENSION(nvar) :: Dif  
!------------------------------ coefficient de diffusion (beta)
      INTEGER :: kpp_ini
!initial condition : 0 => continuous 1 => discontinuous 

!=====================================================

! .....Niveaux max et min, Nbre de points, Nbre iterations
      INTEGER :: niv_max, niv_min, npdt, nbre_iter
! .....Nbre pas de temps, ordre et etendue du support interpolation
      INTEGER :: nt, sgrad, s
! .....Nbre de maille et de points dans l'arbre
      INTEGER :: nbre_maille, Nbre_pt_arbre
! .....Compteurs de feuilles creees, de feuilles reelles et fictives
      INTEGER :: compteur_create, compteur_feuille, compteur_fictive

! .....Nbre de racines des arbres dans les 3 directions
      INTEGER, DIMENSION(3) :: nbre_racine
! .....Limites du domaine et plus petit pas d'espace suivant direction
      DOUBLE PRECISION, DIMENSION(ndim) :: lim_deb, lim_fin, dxmin
! .....Temps, pas de temps, temps et increment des enregistrement, temps fin
      DOUBLE PRECISION :: temps, dt, tps_deb, dtinc, tps_fin
! .....Valeur max du detail, precision, 
      DOUBLE PRECISION :: det_max, epsilon
! .....Min epsilon, max epsilon
      DOUBLE PRECISION :: eps_min, eps_max

! .....Coefficients d'interpolation pour calcul details
      DOUBLE PRECISION, DIMENSION(:), POINTER :: coef

! .... Solution sans MR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: sol_split
      DOUBLE PRECISION, allocatable, dimension(:) :: err_max, error_l2
      INTEGER :: i_err, i_sav, save_res
  
      end module mod_common
