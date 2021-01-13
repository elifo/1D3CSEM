subroutine ReceiverParam (Station,statfile)

use sreceivers

implicit none


type (receiver_stream), intent (INOUT), target :: Station
character (len=20),     intent (IN)            :: statfile


integer :: i
integer, pointer :: ns



open (15, file=statfile, status="old", form="formatted")
read (15,*) Station%nsta

ns => Station%nsta
allocate(Station%Rec(0:ns-1))


! Read coordinates of the station

do  i = 0,ns-1
    read (15,*) Station%Rec(i)%xRec
enddo


! Stress Strain stations

read (15,*)
read (15,*) Station%SSnsta
allocate(Station%SSRec(0:(Station%SSnsta)-1))

do  i = 0,Station%SSnsta-1
    read (15,*) Station%SSRec(i)%SSxRec
enddo
close(15)

return

end subroutine