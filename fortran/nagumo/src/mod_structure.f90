      module mod_structure

      use mod_common

      implicit none

!------------------------------------------------------------

      type structure_fils  
 
        type(structure_maille), pointer :: ptr => null()

      end type structure_fils

!------------------------------------------------------------

      type structure_support
 
        type(structure_maille), pointer :: suivant => null()
        type(structure_maille), pointer :: precedent => null()

      end type structure_support

!------------------------------------------------------------

      type tab_pointer  
 
        double precision, pointer :: t_ptr => null()

      end type tab_pointer



!------------------------------------------------------------	

      type structure_maille 

! .....est-ce que la maille appartient a l'arbre ?
        logical :: arbre

! .....est-ce que la maille est reel (true) ou feuille fictive (false)
        logical :: feuille

! .....niveau de la maille
        integer :: maille_niv 

! .....valeur du detail retenu pour la multiresolution
        real(kind=8) :: detail

! .....indice de la maille dans les directions (ndim)
        integer, dimension(ndim) :: maille_i 

! .....coordonnees et normale de la maille dans les directions (ndim)
        real(kind=8), dimension(ndim) :: x

! .....pas d'espace dans les directions (ndim)
        real(kind=8), dimension(ndim) :: dx

! .....valeur des details par variable (nvar) par rapport aux fils (2**ndim)
        real(kind=8), dimension(nvar, 2**ndim) :: det

! .....flux diffusion a gauche et a droite par variable (nvar)
!      et  par direction (ndim)
        real(kind=8), dimension(nvar, ndim) :: fdiff_gch, fdiff_dt

! .....variables conservatives au temps (n+1)*deltat
        real(kind=8), dimension(nvar) :: u

! .....pointer sur le pere de la maille
        type(structure_maille), pointer :: pere => null()

! .....pointer sur les fils (2**ndim) de la maille
        type(structure_fils), dimension(2**ndim) :: fils

! .....pointer sur les voisins par direction (ndim)
        type(structure_support), dimension(ndim) :: support


      end type structure_maille
!----------------------------------------------------------------
      type structure_feuille 

! .....flux diffusion a gauche et a droite par variable (nvar)
!      et  par direction (ndim)
        type(tab_pointer), dimension(nvar, ndim) :: fdiff_gch, fdiff_dt

! .....pointers pour les variables de la feuille
!        type(tab_pointer) :: var 
        type(tab_pointer), dimension(nvar) :: var 

! .....pas d'espace dans les directions (ndim)
        real(kind=8), dimension(ndim) :: dx

! .....pointer sur les fils (2**ndim) de la maille
        type(structure_fils), dimension(2**ndim) :: fils

! .....pointer sur les voisins par direction (ndim)
        type(structure_support), dimension(ndim) :: support

      end type structure_feuille
!----------------------------------------------------------------



       type(structure_fils), pointer, dimension(:,:,:) :: racine
 

       type(structure_feuille), allocatable, dimension(:) :: tab_phi

      end module mod_structure
