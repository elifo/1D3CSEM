subroutine savefield (Tdomain,GlobCoord,npt,it,TimeP)

use sdomain;use pig
use stimeparam

implicit none

type (domain), intent (IN):: TDomain
type (time), intent (IN), target :: TimeP
integer, intent (IN) :: npt,it
real, dimension  (0:npt-1), intent (IN) :: GlobCoord

! local variables
integer :: ngll,ipoint,i,n,nel
real :: rmodul
real, dimension (0:npt-1) :: Field 
 character(len=20) :: fnamef

nel = Tdomain%nelem

Field = -2e6
do n = 0,nel-1
     ngll = Tdomain%specel(n)%ngllx
     do i = 0,ngll-1
          ipoint = Tdomain%specel(n)%Iglobnum(i)
          Field (ipoint) = Tdomain%specel(n)%Veloc(i,3)
     enddo
enddo
    
! Dump field
write(fnamef,"(a,I6.6)")"velofieldt",it
open (59,file=fnamef,status="unknown",form="formatted")
do ipoint=0,npt-1
  if (.not. Field(ipoint) < -1e6) then
    rmodul=abs(Field(ipoint))
    if (rmodul < 1e-20 ) Field(ipoint) = 0
    write (59,"(4G17.8)") GlobCoord (ipoint), Field(ipoint)
  endif
enddo  
call flush(59); close(59)
return
end subroutine
