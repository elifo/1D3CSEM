subroutine savetrace (Tdomain,Enonli,Station,it,TimeP)

use sdomain
use nonlinear 
use sreceivers
use stimeparam

implicit none 


type (domain),          intent (IN),     target :: Tdomain
type (nonli),           intent (IN),     target :: Enonli
type (receiver_stream), intent (INOUT),  target :: Station
integer,                intent (IN)             :: it
type (time),            intent (IN),     target :: TimeP



integer       :: ir 
integer       :: ni
integer       :: i

real          :: np
real          :: dum0
real          :: dum1
real          :: dum2
real          :: dum3
real          :: dum0x
real          :: dum1x
real          :: dum2x
real          :: dum3x
real          :: dum0y
real          :: dum1y
real          :: dum2y
real          :: dum3y
real          :: dumsat1
real          :: dumsat2
real          :: dumsat3
real          :: dumsat4
real          :: dumsat5
real          :: dumsat6



integer, pointer      :: ngllx
integer, pointer      :: ngllz



if (mod((it-1),TimeP%output) == 0) then
! Output time factor

  ! Acceleration & Velocity histories
  do ir = 0, Station%nsta-1

    ni = Station%Rec(ir)%nr

    dum0 = 0 ; dum0x = 0   ; dum0y = 0
    dum1 = 0 ; dum1x = 0   ; dum1y = 0

    ngllx => Tdomain%specel(ni)%ngllx
    do i =0,ngllx-1 

      dum0 = dum0 + Station%Rec(ir)%Interp_Coeff(i) * Tdomain%specel(ni)%Accel(i,3)    
      dum1 = dum1 + Station%Rec(ir)%Interp_Coeff(i) * Tdomain%specel(ni)%Veloc(i,3)

      if (abs(dum0) < 1e-20) dum0 = 0.0
      if (abs(dum1) < 1e-20) dum1 = 0.0

      dum0x = dum0x + Station%Rec(ir)%Interp_Coeff(i) * Tdomain%specel(ni)%Accel(i,1)    
      dum1x = dum1x + Station%Rec(ir)%Interp_Coeff(i) * Tdomain%specel(ni)%Veloc(i,1)
   
      if (abs(dum0x) < 1e-20) dum0x = 0.0
      if (abs(dum1x) < 1e-20) dum1x = 0.0

      dum0y = dum0y + Station%Rec(ir)%Interp_Coeff(i) * Tdomain%specel(ni)%Accel(i,2)    
      dum1y = dum1y + Station%Rec(ir)%Interp_Coeff(i) * Tdomain%specel(ni)%Veloc(i,2)
   
      if (abs(dum0y) < 1e-20) dum0y = 0.0
      if (abs(dum1y) < 1e-20) dum1y = 0.0
    enddo

    Station%StoreTraceAzz(ir,(it-1)/TimeP%output) = dum0
    Station%StoreTraceVzz(ir,(it-1)/TimeP%output) = dum1

    Station%StoreTraceAxz(ir,(it-1)/TimeP%output) = dum0x
    Station%StoreTraceVxz(ir,(it-1)/TimeP%output) = dum1x

    Station%StoreTraceAyz(ir,(it-1)/TimeP%output) = dum0y
    Station%StoreTraceVyz(ir,(it-1)/TimeP%output) = dum1y
  enddo
endif



if (mod((it-1),TimeP%output) == 0) then

  ! Stress-strain recordings
  do ir = 0, Station%SSnsta-1
    
    ni = Station%SSRec(ir)%SSnr
    
    dum2 = 0  ; dum2x = 0  ; dum2y = 0
    dum3 = 0  ; dum3x = 0  ; dum3y = 0
    
    ngllx => Tdomain%specel(ni)%ngllx
    do i =0,ngllx-1
      dum2 = dum2 + Station%SSRec(ir)%SSInterp_Coeff(i) * Tdomain%specel(ni)%gamma2(i,6)   
      dum3 = dum3 + Station%SSRec(ir)%SSInterp_Coeff(i) * Tdomain%specel(ni)%Stress(i,6)
    
      if (abs(dum2) < 1e-20) dum2 = 0.0
      if (abs(dum3) < 1e-20) dum3 = 0.0

      dum2x = dum2x + Station%SSRec(ir)%SSInterp_Coeff(i) * Tdomain%specel(ni)%gamma2(i,4)   
      dum3x = dum3x + Station%SSRec(ir)%SSInterp_Coeff(i) * Tdomain%specel(ni)%Stress(i,4)

     
      if (abs(dum2x) < 1e-20) dum2x = 0.0
      if (abs(dum3x) < 1e-20) dum3x = 0.0

      dum2y = dum2y + Station%SSRec(ir)%SSInterp_Coeff(i) * Tdomain%specel(ni)%gamma2(i,5)   
      dum3y = dum3y + Station%SSRec(ir)%SSInterp_Coeff(i) * Tdomain%specel(ni)%Stress(i,5)
      
      if (abs(dum2y) < 1e-20) dum2y = 0.0
      if (abs(dum3y) < 1e-20) dum3y = 0.0
    enddo

    Station%SSStoreTrace1zz(ir,(it-1)/TimeP%output) = dum2
    Station%SSStoreTrace2zz(ir,(it-1)/TimeP%output) = dum3

    Station%SSStoreTrace1xz(ir,(it-1)/TimeP%output) = dum2x
    Station%SSStoreTrace2xz(ir,(it-1)/TimeP%output) = dum3x

    Station%SSStoreTrace1yz(ir,(it-1)/TimeP%output) = dum2y
    Station%SSStoreTrace2yz(ir,(it-1)/TimeP%output) = dum3y
  enddo
endif



if (mod((it-1),TimeP%output) == 0) then

  ! Moduli and Poisson ratio
  do ir = 0, Station%SSnsta-1

    ni = Station%SSRec(ir)%SSnr
    
    dumsat1 = 0.0
    dumsat2 = 0.0
    dumsat3 = 0.0
    dumsat4 = 0.0
    dumsat5 = 0.0
    dumsat6 = 0.0


    ngllx => Tdomain%specel(ni)%ngllx
    do i = 0,ngllx-1
      dumsat1 = dumsat1 + Station%SSRec(ir)%SSInterp_Coeff(i) * Tdomain%specel(ni)%Gmodulus  
      dumsat2 = dumsat2 + Station%SSRec(ir)%SSInterp_Coeff(i) * Tdomain%specel(ni)%Emodulus  
      dumsat3 = dumsat3 + Station%SSRec(ir)%SSInterp_Coeff(i) * Tdomain%specel(ni)%Kmodulus
      dumsat4 = dumsat4 + Station%SSRec(ir)%SSInterp_Coeff(i) * Tdomain%specel(ni)%Ni
      dumsat5 = dumsat5 + Station%SSRec(ir)%SSInterp_Coeff(i) * Enonli%specel(ni)%Gm0(i)
      dumsat6 = dumsat6 + Station%SSRec(ir)%SSInterp_Coeff(i) * Enonli%specel(ni)%Gmod2(i)


    
      if (abs(dumsat1) < 1e-20) dumsat1 = 0.0
      if (abs(dumsat2) < 1e-20) dumsat2 = 0.0
      if (abs(dumsat3) < 1e-20) dumsat3 = 0.0
      if (abs(dumsat4) < 1e-20) dumsat4 = 0.0
      if (abs(dumsat5) < 1e-20) dumsat5 = 0.0
      if (abs(dumsat6) < 1e-20) dumsat6 = 0.0

    enddo

    Station%StoreGmod (ir,(it-1)/TimeP%output) = dumsat1
    Station%StoreEmod (ir,(it-1)/TimeP%output) = dumsat2
    Station%StoreKmod (ir,(it-1)/TimeP%output) = dumsat3
    Station%StoreNi   (ir,(it-1)/TimeP%output) = dumsat4
    Station%StoreGm0  (ir,(it-1)/TimeP%output) = dumsat5
    Station%StoreGact (ir,(it-1)/TimeP%output) = dumsat6

  enddo
endif




if (mod((it-1),TimeP%output) == 0  .AND.  Tdomain%Epressmod) then

  ! Saturation Parameters I
  do ir = 0, Station%SSnsta-1

    ni = Station%SSRec(ir)%SSnr
    
    dumsat1 = 0.0
    dumsat2 = 0.0
    dumsat3 = 0.0
    dumsat4 = 0.0
    dumsat5 = 0.0
    dumsat6 = 0.0


    ngllx => Tdomain%specel(ni)%ngllx
    do i = 0,ngllx-1
      dumsat1 = dumsat1 + Station%SSRec(ir)%SSInterp_Coeff(i) * Enonli%specel(ni)%qStress2(i) 
      dumsat2 = dumsat2 + Station%SSRec(ir)%SSInterp_Coeff(i) * Enonli%specel(ni)%r1(i)  
      dumsat3 = dumsat3 + Station%SSRec(ir)%SSInterp_Coeff(i) * Enonli%specel(ni)%Ws2(i)
      dumsat4 = dumsat4 + Station%SSRec(ir)%SSInterp_Coeff(i) * Enonli%specel(ni)%Wn(i)
      dumsat5 = dumsat5 + Station%SSRec(ir)%SSInterp_Coeff(i) * Enonli%specel(ni)%S0(i)
      dumsat6 = dumsat6 + Station%SSRec(ir)%SSInterp_Coeff(i) * Enonli%specel(ni)%S(i)

    
      if (abs(dumsat1) < 1e-20) dumsat1 = 0.0
      if (abs(dumsat2) < 1e-20) dumsat2 = 0.0
      if (abs(dumsat3) < 1e-20) dumsat3 = 0.0
      if (abs(dumsat4) < 1e-20) dumsat4 = 0.0
      if (abs(dumsat5) < 1e-20) dumsat5 = 0.0
      if (abs(dumsat6) < 1e-20) dumsat6 = 0.0
    enddo

    Station%StoreT   (ir,(it-1)/TimeP%output) = dumsat1
    Station%Storer   (ir,(it-1)/TimeP%output) = dumsat2
    Station%StoreWs  (ir,(it-1)/TimeP%output) = dumsat3
    Station%StoreWn  (ir,(it-1)/TimeP%output) = dumsat4
    Station%StoreS0  (ir,(it-1)/TimeP%output) = dumsat5
    Station%StoreS   (ir,(it-1)/TimeP%output) = dumsat6
  enddo
endif



if (mod((it-1),TimeP%output) == 0  .AND.  Tdomain%Epressmod) then

  ! Saturation Parameters II
  do ir = 0, Station%SSnsta-1

    ni = Station%SSRec(ir)%SSnr
    
    dumsat1 = 0.0
    dumsat2 = 0.0
    dumsat3 = 0.0
    dumsat4 = 0.0
    dumsat5 = 0.0


    ngllx => Tdomain%specel(ni)%ngllx
    do i = 0,ngllx-1
      dumsat1 = dumsat1 + Station%SSRec(ir)%SSInterp_Coeff(i) * Enonli%specel(ni)%Pore2(i) 
      dumsat5 = dumsat5 + Station%SSRec(ir)%SSInterp_Coeff(i) * Enonli%specel(ni)%gammaref(i)

    
      if (abs(dumsat1) < 1e-20) dumsat1 = 0.0
      if (abs(dumsat2) < 1e-20) dumsat2 = 0.0
      if (abs(dumsat3) < 1e-20) dumsat3 = 0.0
      if (abs(dumsat4) < 1e-20) dumsat4 = 0.0
      if (abs(dumsat5) < 1e-20) dumsat5 = 0.0
    enddo

    Station%StorePore       (ir,(it-1)/TimeP%output) = dumsat1
    Station%StoreGref       (ir,(it-1)/TimeP%output) = dumsat5

  enddo
endif

return

end subroutine
