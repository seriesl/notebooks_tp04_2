      subroutine flux_euler(dir, ifirst, ilast, u_n, feuler )

      use mod_common

      IMPLICIT NONE

! +++++Variables arguments

      INTEGER, INTENT(IN) :: dir
      INTEGER, INTENT(IN) :: ifirst, ilast

      REAL(kind=8), INTENT(IN), DIMENSION(nvar, ifirst:ilast) :: u_n
      REAL(kind=8), INTENT(OUT), DIMENSION(nvar, ifirst:ilast) ::  feuler

! +++++Variables locales

      INTEGER :: l, m

! ======================= Flux Euler =============================

! .....Cas lineaire : convection a vitesse constante = 1.
      do l = ifirst, ilast
        do m = 1, nvar
          feuler(m,l) = 0.5*u_n(dir,l)**2
        enddo
      enddo

      end subroutine flux_euler
