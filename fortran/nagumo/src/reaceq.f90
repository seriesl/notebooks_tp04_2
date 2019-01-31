      subroutine reaceq (ni,t,pphi,dphidt,rpar,ipar)

      use mod_common

      implicit none

      integer, intent(in) ::  ni, ipar
      double precision, intent(in)  ::  t, rpar
      double precision,dimension(ni), intent(in)   :: pphi
      double precision,dimension(ni), intent(out)  :: dphidt

      double precision,dimension(ni)   :: phi      

         if (pphi(1) < 0.d0) then
             phi(1) = 0.d0
         else
             phi(1) = pphi(1)
         endif
 
      dphidt(1) =  (1.d0/Dif(1))*phi(1)*phi(1)*(1.d0 - phi(1))

      end subroutine reaceq

