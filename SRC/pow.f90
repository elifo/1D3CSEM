real function pow (x,vp,npm,dx,np,A)

implicit none

real :: x,vp,dx,A
integer :: npm,np

! Local variables

real :: rnpm,pp1

rnpm = float (npm)


pp1 = dx
pp1 = 1 /pp1
pow = A * vp * pp1 * (x/rnpm)**np

return
end function
