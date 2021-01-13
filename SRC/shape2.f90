subroutine shape2(Tdomain,nm,GlobCoord)

use sdomain

implicit none

type(domain),target, intent (INOUT) :: Tdomain
integer, intent (IN) :: nm
real, intent (INOUT) :: GlobCoord(0:nm)

integer :: nel,i,n,np,i_domain
integer, pointer :: ngllx
real :: xiv
real, dimension (0:1) :: vcoord

nel = Tdomain%nelem

do n = 0,nel - 1
   vcoord(0) = Tdomain%node(Tdomain%specel(n)%Enode(0)-1)%coord
   vcoord(1) = Tdomain%node(Tdomain%specel(n)%Enode(1)-1)%coord
   ngllx => Tdomain%specel(n)%ngllx
   i_domain = Tdomain%specel(n)%Numdomain-1
   do i = 0,ngllx - 1
     xiv = Tdomain%Sub_Domain(i_domain)%GLLc(i) 
     np = Tdomain%specel(n)%Iglobnum(i)
     GlobCoord(np) = 0.5*(vcoord(0)*(1-xiv)+vcoord(1)*(1+xiv))
   enddo 
         
enddo

return
end subroutine
