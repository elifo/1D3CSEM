subroutine velocdiric (vmax,t0,tmax,timed,vb)

use sdomain
use stimeparam

implicit none

real, intent (IN) :: timed,vmax,t0,tmax
real, intent (INOUT) :: vb


if ( timed < t0 )  then
   vb = (vmax/t0) * timed
   vb = 0.
endif
if ( timed >= t0 .and. timed <= (tmax-t0) ) then
   vb = vmax
endif
if ( timed > (tmax-t0) ) then
   vb = (vmax/t0) * (-timed+tmax)
   vb = 0.
endif

!if (vmax > 0) write(25,*) timed,vb


return
end subroutine velocdiric
