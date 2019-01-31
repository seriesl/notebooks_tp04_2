      subroutine diffeq (ni,t,phi,dphidt,rpar,ipar)

      use mod_common
      use mod_structure

      implicit none

      integer, intent(in) ::  ni, ipar
      double precision, intent(in)  ::  t, rpar
      double precision,dimension(ni), intent(inout)   :: phi
      double precision,dimension(ni), intent(out)  :: dphidt
      
      integer ::  i


      call fonction_diffusion(ni, phi, dphidt, rpar, ipar)



      end subroutine diffeq

