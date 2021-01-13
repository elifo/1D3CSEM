subroutine SourceParameters(Source,sourcefile)

use ssource

implicit none

type(sour), target, intent (INOUT) :: Source
character(len=20),  intent (IN)    :: sourcefile


! Read coordinates of the source

open (14, file=sourcefile, status="old", form="formatted")
read (14,*) Source%xsp


! Read type of Source : the type is an integer with the following
! meaning : 1 Impulse Source 
!           2 Explosion

read (14,*) Source%kind


! Read function form of the Source: 1 Gaussian
!                                   2 Ricker

read (14,*) Source%funct



if (Source%funct == 1) then
	read (14,*) Source%tau    ! Width for Gaussian
    read (14,*) Source%f0
else if (Source%funct == 2) then
    read (14,*) Source%tau    ! Delay time for Ricker wave
    read (14,*) Source%f0     ! Cut off frequency

else
    write(*,"(A)", ADVANCE = "NO") "Not supported functional form for the source - Press <Return> to continue."
    read(*,"(A)") 
endif


return
end subroutine 