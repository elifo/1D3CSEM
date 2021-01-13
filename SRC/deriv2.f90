subroutine deriv2 (Tdomain)

use sdomain

implicit none

type (domain), intent (INOUT), target :: Tdomain

integer :: nel,i,n
integer, pointer :: ngllx
real :: Jac
real, dimension (0:1) :: vcoord

! This subroutine calculates the derivatives of the shape function and
! defines the DirectGrad 

nel  = Tdomain%nelem

do n = 0 ,nel - 1

   vcoord(0) = Tdomain%node(Tdomain%specel(n)%Enode(0)-1)%coord
   vcoord(1) = Tdomain%node(Tdomain%specel(n)%Enode(1)-1)%coord
   ngllx => Tdomain%specel(n)%ngllx
   Jac = 0.5*(vcoord(1)-vcoord(0))
   do i = 0, ngllx - 1       
      Tdomain%specel(n)%Jacob(i) = Jac
      Tdomain%specel(n)%InvGrad(i) = 1/Jac
   enddo     
enddo
return
end subroutine
