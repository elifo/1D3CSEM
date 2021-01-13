subroutine ReceiverPosition(Tdomain,Station,GlobCoord,npt)

use sdomain
use sreceivers

implicit none  

type (receiver_stream), intent (INOUT), target :: Station
type (domain),          intent (IN),    target :: Tdomain
integer,                intent (IN)            :: npt
real, dimension (0:npt-1) :: GlobCoord


integer      :: n
integer      :: nrec
integer      :: i 
integer      :: nwh
integer      :: nglxinput
integer      :: i_domain
integer      :: node0
integer      :: node1
integer      :: n_element

real         :: x0 
real         :: x1 
real         :: outx
real         :: xi1
real         :: dwdxi
real         :: wxi

integer, pointer :: ngllx
integer, pointer :: nel 
integer, pointer :: ns 


nel   => Tdomain%nelem
ns    => Station%nsta



! Compute the nearest point to the real location of the receivers in GLL scheme

do nrec = 0,ns-1

  ! Search for the element
  do_rec1 : do n = 0,nel-1
    node0 = Tdomain%specel(n)%Enode(0)-1
    node1 = Tdomain%specel(n)%Enode(1)-1
    x0 = Tdomain%node(node0)%coord 
    x1 = Tdomain%node(node1)%coord 
    if (x0 <= Station%Rec(nrec)%xRec .and. Station%Rec(nrec)%xRec <= x1 ) then
      n_element = n
      exit do_rec1
    endif
  enddo do_rec1
 

  ! Determine coordinates 
  Station%Rec(nrec)%nr = n_element 

  x0  = Tdomain%node(Tdomain%specel(n_element)%Enode(0)-1)%coord
  x1  = Tdomain%node(Tdomain%specel(n_element)%Enode(1)-1)%coord
  xi1 = ( 2*Station%Rec(nrec)%xRec - (x0+x1) ) / (x1-x0)  

  Station%Rec(nrec)%xi = xi1
  i_domain  = Tdomain%specel(n)%Numdomain-1
  nglxinput = Tdomain%specel(Station%Rec(nrec)%nr)%ngllx
  allocate(Station%Rec(nrec)%Interp_Coeff(0:nglxinput-1))

  do i = 0,nglxinput -1
    call  pol_lagrange(nglxinput,Tdomain%Sub_Domain(i_domain)%GLLc,i,Station%Rec(nrec)%xi,outx)
    Station%Rec(nrec)%Interp_Coeff(i) = outx
  enddo
enddo 

!
! Stress- strain receiver interpolation 
!
ns => Station%SSnsta
! Compute the nearest point to the real location of the receivers in GLL scheme

do nrec = 0,ns-1
  ! Search for the element
  do_SSrec1 : do n = 0,nel-1
    node0 = Tdomain%specel(n)%Enode(0)-1
    node1 = Tdomain%specel(n)%Enode(1)-1

    x0 = Tdomain%node(node0)%coord 
    x1 = Tdomain%node(node1)%coord

    if (x0 <= Station%SSRec(nrec)%SSxRec .and. Station%SSRec(nrec)%SSxRec <= x1) then
      n_element = n
      exit do_SSrec1
    endif
  enddo do_SSrec1
 
  ! Determine coordinates 
  Station%SSRec(nrec)%SSnr = n_element 

  x0  = Tdomain%node(Tdomain%specel(n_element)%Enode(0)-1)%coord
  x1  = Tdomain%node(Tdomain%specel(n_element)%Enode(1)-1)%coord
  xi1 = ( 2*Station%SSRec(nrec)%SSxRec - (x0+x1) ) / (x1-x0)  
  Station%SSRec(nrec)%SSxi = xi1

  i_domain  = Tdomain%specel(n)%Numdomain-1
  nglxinput = Tdomain%specel(Station%SSRec(nrec)%SSnr)%ngllx
  allocate(Station%SSRec(nrec)%SSInterp_Coeff(0:nglxinput-1))


  do i = 0,nglxinput -1
    call  pol_lagrange(nglxinput,Tdomain%Sub_Domain(i_domain)%GLLc,i,Station%SSRec(nrec)%SSxi,outx)
    Station%SSRec(nrec)%SSInterp_Coeff(i) = outx
  enddo


enddo 








return

end subroutine