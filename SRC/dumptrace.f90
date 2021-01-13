subroutine dumptrace (Station,TimeP,Tdomain)

use sreceivers
use stimeparam
use sdomain


implicit none

type (receiver_stream), intent (IN), target :: Station
type (time),            intent (IN), target :: TimeP
type (domain),          intent (IN), target :: Tdomain



integer     :: i
integer     :: i1
integer     :: it 
integer     :: nt 

real        :: t
real        :: dum
real        :: max1
real        :: max2 
real        :: deltat

character*40   :: fnamef




! Dumping the traces

nt      = TimeP%nt
deltat  = TimeP%dt


do i = 0,Station%nsta-1
  i1 = i+1
  
  ! xz direction
  write (fnamef,"(a,I3.3)") 'outputfiles/accelx',i1
  open (31,file=fnamef, form="formatted")

  write (fnamef,"(a,I3.3)") 'outputfiles/velocx',i1
  open (32,file=fnamef, form="formatted")

  do it = 1, nt
    t = float(it-1)*deltat
    
    if (mod((it-1),TimeP%output)==0) then
      write (31,*) t, Station%StoreTraceAxz (i,(it-1)/TimeP%output)
      write (32,*) t, Station%StoreTraceVxz (i,(it-1)/TimeP%output)
    endif
  enddo
  close (31); close(32)

  ! yz direction
  write (fnamef,"(a,I3.3)") 'outputfiles/accely',i1
  open (33,file=fnamef, form="formatted")

  write (fnamef,"(a,I3.3)") 'outputfiles/velocy',i1
  open (34,file=fnamef, form="formatted")

  do it = 1, nt
    t = float(it-1)*deltat
    
    if (mod((it-1),TimeP%output)==0) then
      write (33,*) t, Station%StoreTraceAyz (i,(it-1)/TimeP%output)
      write (34,*) t, Station%StoreTraceVyz (i,(it-1)/TimeP%output)
    endif
  enddo
  close (33); close(34)

  ! zz direction
  write (fnamef,"(a,I3.3)") 'outputfiles/accelz',i1
  open (35,file=fnamef, form="formatted")

  write (fnamef,"(a,I3.3)") 'outputfiles/velocz',i1
  open (36,file=fnamef, form="formatted")

  do it = 1, nt
    t = float(it-1)*deltat
    
    if (mod((it-1),TimeP%output)==0) then
      write (35,*) t, Station%StoreTraceAzz (i,(it-1)/TimeP%output)
      write (36,*) t, Station%StoreTraceVzz (i,(it-1)/TimeP%output)
    endif
  enddo
  close (35); close(36)
enddo



do i = 0,Station%SSnsta-1
  i1 = i+1

  ! Stress strain 
  ! xz direction
  write (fnamef,"(a,I3.3)") 'outputfiles/StressStrainxz',i1
  open (61,file=fnamef, form="formatted")

  do it = 1, nt
    t = float(it-1)*deltat

    if (mod((it-1),TimeP%output)==0) then
      write (61,*) Station%SSStoreTrace1xz(i,(it-1)/TimeP%output), &
                   Station%SSStoreTrace2xz(i,(it-1)/TimeP%output)
    endif
  enddo
  close (61)

  ! yz direction
  write (fnamef,"(a,I3.3)") 'outputfiles/StressStrainyz',i1
  open (62,file=fnamef, form="formatted")

  do it = 1, nt
    t = float(it-1)*deltat

    if (mod((it-1),TimeP%output)==0) then
      write (62,*) Station%SSStoreTrace1yz(i,(it-1)/TimeP%output), &
                   Station%SSStoreTrace2yz(i,(it-1)/TimeP%output)
    endif
  enddo
  close (62)

  ! zz direction
  write (fnamef,"(a,I3.3)") 'outputfiles/StressStrainzz',i1
  open (63,file=fnamef, form="formatted")

  do it = 1, nt
    t = float(it-1)*deltat

    if (mod((it-1),TimeP%output)==0) then
      write (63,*) Station%SSStoreTrace1zz(i,(it-1)/TimeP%output), &
                   Station%SSStoreTrace2zz(i,(it-1)/TimeP%output)
    endif
  enddo
  close (63)



if (Tdomain%Epressmod) then
  ! Soil Parameters
  write (fnamef,"(a,I3.3)") 'outputfiles/PressSoilParams',i1
  open (64,file=fnamef, form="formatted")

  do it = 1, nt
    t = float(it-1)*deltat

    if (mod((it-1),TimeP%output)==0) then
      write (64,*) Station%StoreGmod(i,(it-1)/TimeP%output), &
                   Station%StoreEmod(i,(it-1)/TimeP%output), &
                   Station%StoreKmod(i,(it-1)/TimeP%output), &
                   Station%StoreNi  (i,(it-1)/TimeP%output), &
                   Station%StoreGm0 (i,(it-1)/TimeP%output), &
                   Station%StoreGact(i,(it-1)/TimeP%output)

    endif
  enddo
  close (64)
endif



if (Tdomain%Epressmod) then
  ! Effective Stress Analysis Parameters 
  write (fnamef,"(a,I3.3)") 'outputfiles/PressEffectiveParams',i1
  open (65,file=fnamef, form="formatted")

  do it = 1, nt
    t = float(it-1)*deltat

    if (mod((it-1),TimeP%output)==0) then
      write (65,*) Station%StoreT     (i,(it-1)/TimeP%output), &
                   Station%Storer     (i,(it-1)/TimeP%output), &
                   Station%StoreWs    (i,(it-1)/TimeP%output), &
                   Station%StoreWn    (i,(it-1)/TimeP%output), &
                   Station%StoreS0    (i,(it-1)/TimeP%output), &
                   Station%StoreS     (i,(it-1)/TimeP%output), &
                   Station%StorePore  (i,(it-1)/TimeP%output), &
                   Station%StoreGref  (i,(it-1)/TimeP%output)
    endif
  enddo
  close (65)
endif



enddo

return

end subroutine dumptrace