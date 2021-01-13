subroutine TestCFL(Tdomain,TimeP,nm,GlobCoord)

use sdomain
use stimeparam


implicit none


type (domain), intent(IN),  target :: Tdomain
type (time),   intent (IN), target :: TimeP
integer,       intent (IN)         :: nm
real,          intent (IN)         :: GlobCoord(0:nm)



integer        :: n
integer        :: np1
integer        :: np2
integer        :: i_domain

real           :: vp
real           :: dx
real           :: cfl
real           :: cflmax
real           :: cfllimit
real           :: dxmin

integer, pointer  :: nel




nel => Tdomain%nelem

cflmax   = 0.0 
dxmin    = 10000.0
cfllimit = 0.3


do n = 0,nel-1

    ! Minimum distance between GLL points
    np1 = Tdomain%specel(n)%Iglobnum(1)
    np2 = Tdomain%specel(n)%Iglobnum(0)
    dx  = abs(GlobCoord(np1) - GlobCoord(np2)) 

    i_domain = Tdomain%specel(n)%NumDomain-1
!   vp  = Tdomain%Sub_domain(i_domain)%Mu/ Tdomain%Sub_domain(i_domain)%rho

    vp  = Tdomain%Sub_domain(i_domain)%lambda2mu/ Tdomain%Sub_domain(i_domain)%rho
    vp  = sqrt(vp)
    cfl = vp*Tdomain%specel(n)%dt/dx

    if (cfl > cflmax)  cflmax = cfl
    if (dx  < dxmin)   dxmin  = dx     
enddo

print*,'CFL   ',  cfl
print*,'dxmin ',  dxmin


if ( cfl > cfllimit ) then
    write(*,"(A)", ADVANCE = "NO") "CFL greater than CFLmax = 0.3 - Press <Return> to continue."
	read (*,"(A)")
endif

return
end subroutine TestCFL