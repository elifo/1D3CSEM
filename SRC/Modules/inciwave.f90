module inciwave

use astdio

type :: inci
 real, pointer, dimension(:,:) :: VB			!!!
 character (len = 20) :: imposig
 logical :: imposource
 integer :: tfactor, sposi
 real    :: timeup, coeff
end type


contains


subroutine incidentw(Ewave,dtsim)
! SOURCE IMPOSITION by VELOCITY 
!

implicit none

type(inci),intent(INOUT),target :: Ewave
real, intent(IN), target :: dtsim

real :: elo,dtwave
integer :: j, N
real, dimension(:), allocatable :: nt

print*, "Incident wave velocity assignment"
N = IO_file_length(Ewave%imposig)


allocate( Ewave%VB(N,3) )
allocate( nt(N))
Ewave%VB = 0.0

open(18,file=Ewave%imposig,status='old')
do j= 1,N
	read (18,*) nt(j), Ewave%VB(j,1), Ewave%VB(j,2), Ewave%VB(j,3)
end do
close(18) 

Ewave%timeup = nt(N)

dtwave = nt(2)-nt(1) 
! print*, dtwave, dtsim, dtwave/dtsim
Ewave%tfactor = int(dtwave* 1.00001/ dtsim)


print*, "Time step factor of input and simulation"
print*, Ewave%tfactor




end subroutine incidentw

end module inciwave
