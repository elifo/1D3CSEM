subroutine SourcePosition(Tdomain,Source,GlobCoord,npt)

use sdomain; 
use ssource; 

implicit none  

type (sour),   intent (INOUT), target :: Source
type (domain), intent (IN),    target :: Tdomain
integer,       intent (IN)            :: npt


integer :: n
integer :: i
integer :: nwh
integer :: nind
integer :: i_domain
integer :: node1
integer :: node0

real    :: x0
real    :: x1 
real    :: xi1
real    :: dwdxi
real    :: wxi

integer, pointer :: ngllx
integer, pointer :: nel 

integer, dimension (0:1) :: nsource

real,    dimension (0:npt-1) :: GlobCoord

real,    dimension (:), pointer :: Jac
real,    dimension (:), pointer :: LocInvGrad




nel => Tdomain%nelem

! Compute the nearest point to the real location of the source in GLL scheme

! Search for the element
nind = 0
do n = 0,nel-1
  ngllx => Tdomain%specel(n)%ngllx

  node0 = Tdomain%specel(n)%Enode(0)-1
  node1 = Tdomain%specel(n)%Enode(1)-1

  x0 = Tdomain%node(node0)%coord 
  x1 = Tdomain%node(node1)%coord 

  if (x0 <= Source%xsp .and. Source%xsp <= x1 ) then
    nsource(nind) = n
    nind = nind + 1
    if (nind > 2 ) then
      write(*,"(A)", ADVANCE = "NO") "Maximum allowed connection is 2 - Press <Return> to continue."
      read(*,"(A)")
    endif
  endif

enddo


! Determine coordinates 
Source%ine = nind

allocate(Source%Elem(0:nind-1))
do n = 0,nind-1
  x0  = Tdomain%node(Tdomain%specel(nsource(n))%Enode(0)-1)%coord
  x1  = Tdomain%node(Tdomain%specel(nsource(n))%Enode(1)-1)%coord
  xi1 = ( 2*Source%xsp - (x0+x1) ) / (x1-x0)  
  Source%Elem(n)%nr = nsource (n)
  if (xi1 >  1.0  .and. xi1 <  1+ 1e-5)  xi1 =  1
  if (xi1 < -1.0  .and. xi1 > -1- 1e-5)  xi1 = -1
  Source%Elem(n)%xi = xi1
enddo


if (Source%kind == 1) then        ! Pulse dire 
  do n = 0, Source%ine-1

    nwh   = Source%Elem(n)%nr
    ngllx => Tdomain%specel(nwh)%ngllx
    Jac   => Tdomain%specel(nwh)%Jacob

    allocate(Source%Elem(n)%Impulse(0:ngllx-1))
    i_domain = Tdomain%specel(nwh)%NumDomain-1

    do i = 0,ngllx-1
      call pol_lagrange(ngllx, Tdomain%Sub_Domain(i_domain)%GLLc, i, Source%Elem(n)%xi, wxi ) ! h_i(x_s)
      Source%Elem(n)%Impulse(i) = wxi  
    enddo
  enddo 
else if (Source%kind == 2) then   ! Explosive source diagonal moment considered
  do n = 0, Source%ine-1

    nwh        = Source%Elem(n)%nr
    ngllx      => Tdomain%specel(nwh)%ngllx
    LocInvGrad => Tdomain%specel(nwh)%InvGrad
    Jac        => Tdomain%specel(nwh)%Jacob

    allocate(Source%Elem(n)%Explosion(0:ngllx-1))
    i_domain = Tdomain%specel(nwh)%NumDomain-1

    do i = 0,ngllx-1
      call DERIVAL(Tdomain%Sub_Domain(i_domain)%GLLc, ngllx, i, Source%Elem(n)%xi, dwdxi)
      Source%Elem(n)%Explosion(i) = dwdxi*LocInvGrad(i) 
    enddo
  enddo

endif
return

end subroutine   
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine DERIVAL (GLLc,n,ip,x0,dy0) 

implicit none

integer, intent (IN) :: n,ip
real,    intent (IN) :: x0
real,    intent (OUT):: dy0
real,    dimension (0:n-1), intent (IN) :: GLLc

integer  :: i
integer  :: k
real     :: y0

dy0  = 0
do i = 0,n-1
  y0 = 1
  if (i /= ip) then
    do k = 0,n-1
      if (k /= i .and. k /= ip)  then
        y0  = y0 * (x0 - GLLc(k) ) /(GLLc(ip)-GLLc(k) ) 
      endif
    enddo
    y0  = y0 / (GLLc(ip)-GLLc(i) ) 
    dy0 = dy0 + y0
  endif
enddo

return

end subroutine DERIVAL
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc