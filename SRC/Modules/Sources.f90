module selem_source
type :: elem_source
integer :: nr   ! Number of the element
real :: xi      ! Coord in ref square
real, dimension (:), pointer :: Impulse
real, dimension (:), pointer :: Explosion
end type 
end module selem_source

module ssource
use selem_source
type :: sour
real :: xsp       ! Source coord
real :: tau,f0    ! Time function parameters
integer :: kind   ! Impulsion or Explosion       !, direction, ine
integer :: funct  ! Gauss or Ricker
integer :: ine    ! Number of elements
type (elem_source), dimension (:),pointer  :: Elem

end type 

contains

real function CompSource (Source,time)

type (sour) :: Source
real :: time

!select case (Source%kind)
!case (1)   ! Impulse
     select case (Source%funct)
        case (1) 
                CompSource = Gaussian (time,Source%tau,Source%f0)
        case (2)
                CompSource = Ricker (time,Source%tau,Source%f0)
         end select
!      
!case (2) ! Explosion
!        select case (Source%funct)
!	case (1)
!		CompSource = Gaussian (time,Source%tau)
!	case (2)
!		CompSource = Ricker (time,Source%tau,Source%f0)
!
!        end select
!end select
return
end function

real function Gaussian (time, tau,f0)

real :: tau,time
Gaussian = -2*(time-tau)*f0**2 * exp (-(time-tau)**2*f0**2)
!Gaussian = exp (-(time-tau)**2*f0**2)  !-(time-tau) * exp (-(time-tau)**2/tau**2)
write(22,*) time,Gaussian

return
end function

real function Ricker (time,tau,f0)

use pig

real :: time, tau, f0

real :: sigma

sigma = pi * f0 * (time - tau )
sigma = sigma **2

Ricker = (1-2 * sigma) *exp (-sigma)
write(23,*) time,Ricker

return
end function



end module
