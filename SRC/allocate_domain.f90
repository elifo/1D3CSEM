subroutine allocate_domain (Tdomain,Evisco,Enonli,Station,eldim)

use sdomain
use logical_input 
use viscopara
use nonlinear 
use sreceivers


implicit none

type(domain), intent (INOUT), target :: Tdomain
type(visco) , intent (INOUT), target :: Evisco
type(nonli) , intent (INOUT), target :: Enonli
type(receiver_stream), intent (INOUT), target :: Station
integer,      intent(IN) :: eldim


integer  :: n
integer, pointer :: nel
integer, pointer :: ngllx


! Receiver allocation

allocate(Station%StoreTraceAzz      (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreTraceAxz      (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreTraceAyz      (0:Station%nsta- 1,0:eldim- 1))

allocate(Station%StoreTraceVzz      (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreTraceVyz      (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreTraceVxz      (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%SSStoreTrace1zz    (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%SSStoreTrace1xz    (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%SSStoreTrace1yz    (0:Station%nsta- 1,0:eldim- 1))

allocate(Station%SSStoreTrace2zz    (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%SSStoreTrace2xz    (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%SSStoreTrace2yz    (0:Station%nsta- 1,0:eldim- 1))

allocate(Station%StoreGmod    (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreEmod    (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreKmod    (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreNi      (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreGm0     (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreT       (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%Storer       (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreWs      (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreWn      (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreS0      (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreS       (0:Station%nsta- 1,0:eldim- 1))

allocate(Station%StorePore    (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreSOcto   (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreGOcto   (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreGref    (0:Station%nsta- 1,0:eldim- 1))
allocate(Station%StoreGact    (0:Station%nsta- 1,0:eldim- 1))



nel => Tdomain%nelem
allocate(Enonli%specel(0:nel-1))


do n = 0,nel-1
    ngllx => Tdomain%specel(n)%ngllx
    allocate(Tdomain%specel(n)%MassMat(0:ngllx-1))
  
    ! Allocation for all the mesh               
    allocate(Tdomain%specel(n)%Stress   (0:ngllx-1,6))    
    allocate(Tdomain%specel(n)%Veloc    (0:ngllx-1,3))
    allocate(Tdomain%specel(n)%Forces   (0:ngllx-1,3))
    allocate(Tdomain%specel(n)%Acoeff   (0:ngllx-1))
    allocate(Tdomain%specel(n)%Calc0    (0:ngllx-1,3)) 
    allocate(Tdomain%specel(n)%Accel    (0:ngllx-1,3))
    allocate(Tdomain%specel(n)%V0       (0:ngllx-1,3)) 
    allocate(Tdomain%specel(n)%MidDispl (0:ngllx-1,6))
    allocate(Tdomain%specel(n)%gamma1   (0:ngllx-1,6))
    allocate(Tdomain%specel(n)%gamma2   (0:ngllx-1,6))
    allocate(Tdomain%specel(n)%epsm     (0:ngllx-1))
    allocate(Tdomain%specel(n)%sigm     (0:ngllx-1))
    allocate(Tdomain%specel(n)%depsm    (0:ngllx-1))
    allocate(Tdomain%specel(n)%dsigm    (0:ngllx-1))
    allocate(Tdomain%specel(n)%deps     (0:ngllx-1,6))
    allocate(Tdomain%specel(n)%de       (0:ngllx-1,6))
    allocate(Tdomain%specel(n)%dsigma   (0:ngllx-1,6))
    allocate(Tdomain%specel(n)%dplastic (0:ngllx-1,6))
    allocate(Tdomain%specel(n)%SigmaIni (0:ngllx-1,6))  

    Tdomain%specel(n)%Veloc     = 0.0
    Tdomain%specel(n)%Accel     = 0.0
    Tdomain%specel(n)%Calc0     = 0.0
    Tdomain%specel(n)%Forces    = 0.0
    Tdomain%specel(n)%Stress    = 0.0
    Tdomain%specel(n)%V0        = 0.0
    Tdomain%specel(n)%MidDispl  = 0.0
    Tdomain%specel(n)%gamma1    = 0.0
    Tdomain%specel(n)%gamma2    = 0.0
    Tdomain%specel(n)%epsm      = 0.0
    Tdomain%specel(n)%sigm      = 0.0
    Tdomain%specel(n)%depsm     = 0.0
    Tdomain%specel(n)%dsigm     = 0.0
    Tdomain%specel(n)%deps      = 0.0
    Tdomain%specel(n)%de        = 0.0
    Tdomain%specel(n)%dsigma    = 0.0
    Tdomain%specel(n)%dplastic  = 0.0
    Tdomain%specel(n)%SigmaIni  = 0.0 


    ! PML condition
    allocate(Tdomain%specel(n)%DumpMass  (0:ngllx-1))
    allocate(Tdomain%specel(n)%DumpSx1   (0:ngllx-1)) 
    allocate(Tdomain%specel(n)%DumpSx2   (0:ngllx-1)) 
    allocate(Tdomain%specel(n)%DumpVx1   (0:ngllx-1)) 
    allocate(Tdomain%specel(n)%DumpVx2   (0:ngllx-1))

    allocate(Tdomain%specel(n)%DumpMassE (0:ngllx-1))
    allocate(Tdomain%specel(n)%DumpSx1E  (0:ngllx-1)) 
    allocate(Tdomain%specel(n)%DumpSx2E  (0:ngllx-1)) 
    allocate(Tdomain%specel(n)%DumpVx1E  (0:ngllx-1)) 
    allocate(Tdomain%specel(n)%DumpVx2E  (0:ngllx-1))

    Tdomain%specel(n)%DumpMass  = 0.0
    Tdomain%specel(n)%DumpSx1   = 0.0
    Tdomain%specel(n)%DumpSx2   = 0.0
    Tdomain%specel(n)%DumpVx1   = 0.0
    Tdomain%specel(n)%DumpVx2   = 0.0

    Tdomain%specel(n)%DumpMassE = 0.0
    Tdomain%specel(n)%DumpSx1E  = 0.0
    Tdomain%specel(n)%DumpSx2E  = 0.0
    Tdomain%specel(n)%DumpVx1E  = 0.0
    Tdomain%specel(n)%DumpVx2E  = 0.0


    ! Visco Array & Matrix Allocation
    if (.NOT. Tdomain%specel(n)%PML ) then
        allocate(Evisco%specel(n)%DumpSx3    (0:ngllx-1)) 
        allocate(Evisco%specel(n)%viszipt    (0:ngllx-1,8))
        allocate(Evisco%specel(n)%DumpSx3YZ  (0:ngllx-1)) 
        allocate(Evisco%specel(n)%visziptYZ  (0:ngllx-1,8))
        allocate(Evisco%specel(n)%DumpSx3ZZ  (0:ngllx-1)) 
        allocate(Evisco%specel(n)%DumpSx3VOL  (0:ngllx-1)) 
        allocate(Evisco%specel(n)%visziptZZ  (0:ngllx-1,8))
        allocate(Evisco%specel(n)%visziptVOL (0:ngllx-1,8))  

        Evisco%specel(n)%DumpSx3    = 0.0
        Evisco%specel(n)%viszipt    = 0.0
        Evisco%specel(n)%DumpSx3YZ  = 0.0
        Evisco%specel(n)%visziptYZ  = 0.0
        Evisco%specel(n)%DumpSx3ZZ  = 0.0
        Evisco%specel(n)%DumpSx3VOL = 0.0
        Evisco%specel(n)%visziptZZ  = 0.0
        Evisco%specel(n)%visziptVOL = 0.0
    end if 

  
    ! Incident Wave Modifications
    allocate(Tdomain%specel(n)%elStress  (0:ngllx-1)) 
    allocate(Tdomain%specel(n)%elMidDispl(0:ngllx-1)) 
    
    Tdomain%specel(n)%elStress      = 0.0
    Tdomain%specel(n)%elMidDispl    = 0.0


    ! Nonlinearity (old algorithm)
    allocate(Enonli%specel(n)%alf   (0:Enonli%Nspr-1,0:ngllx-1))
    allocate(Enonli%specel(n)%pbar  (2,0:ngllx-1))
    allocate(Enonli%specel(n)%dTau  (0:ngllx-1))
    
    Enonli%specel(n)%alf            = 0.0
    Enonli%specel(n)%pbar           = 0.0
    Enonli%specel(n)%dTau           = 0.0


    if (Tdomain%specel(n)%FPML) then
        allocate(Tdomain%specel(n)%IV       (0:ngllx-1))
        allocate(Tdomain%specel(n)%IS       (0:ngllx-1))
        allocate(Tdomain%specel(n)%IVE      (0:ngllx-1))
        allocate(Tdomain%specel(n)%ISE      (0:ngllx-1))             
        allocate(Tdomain%specel(n)%IVeloc   (0:ngllx-1))
        allocate(Tdomain%specel(n)%IStress  (0:ngllx-1))

        Tdomain%specel(n)%IVeloc    = 0.0 
        Tdomain%specel(n)%IStress   = 0.0
    endif       


    ! Total and effective stress analyses
    allocate(Enonli%specel(n)%hexaSigma  (0:ngllx-1,6))
    allocate(Enonli%specel(n)%hexaS2     (0:ngllx-1,6))
    allocate(Enonli%specel(n)%hexaS1     (0:ngllx-1,6))
    allocate(Enonli%specel(n)%dhexaSigma (0:ngllx-1,6))
    allocate(Enonli%specel(n)%dhexaS     (0:ngllx-1,6))
    allocate(Enonli%specel(n)%qStress1   (0:ngllx-1))
    allocate(Enonli%specel(n)%qStress2   (0:ngllx-1))
    allocate(Enonli%specel(n)%p0Stress   (0:ngllx-1))
    allocate(Enonli%specel(n)%Gm0        (0:ngllx-1))
    allocate(Enonli%specel(n)%S0         (0:ngllx-1))
    allocate(Enonli%specel(n)%S          (0:ngllx-1))
    allocate(Enonli%specel(n)%gammaref   (0:ngllx-1))
    allocate(Enonli%specel(n)%eliflogic  (0:ngllx-1))
    allocate(Enonli%specel(n)%Wn         (0:ngllx-1))
    allocate(Enonli%specel(n)%r0         (0:ngllx-1))
    allocate(Enonli%specel(n)%r1         (0:ngllx-1))
    allocate(Enonli%specel(n)%S0old      (0:ngllx-1))
    allocate(Enonli%specel(n)%Sold       (0:ngllx-1))    
    allocate(Enonli%specel(n)%Ws1        (0:ngllx-1))
    allocate(Enonli%specel(n)%Ws2        (0:ngllx-1))
    allocate(Enonli%specel(n)%SigmaP     (0:ngllx-1,3)) 
    allocate(Enonli%specel(n)%Gmod1      (0:ngllx-1))
    allocate(Enonli%specel(n)%Gmod2      (0:ngllx-1))
    allocate(Enonli%specel(n)%Emod1      (0:ngllx-1))
    allocate(Enonli%specel(n)%Emod2      (0:ngllx-1))    
    allocate(Enonli%specel(n)%Kmod1      (0:ngllx-1))
    allocate(Enonli%specel(n)%Kmod2      (0:ngllx-1))    
    allocate(Enonli%specel(n)%gammaSAT   (Enonli%Nspr,  0:ngllx-1))
    allocate(Enonli%specel(n)%RSAT       (Enonli%Nspr,  0:ngllx-1))
    allocate(Enonli%specel(n)%CNinvSAT   (Enonli%Nspr-1,0:ngllx-1))    
    allocate(Enonli%specel(n)%SigmaMEff1 (0:ngllx-1))
    allocate(Enonli%specel(n)%SigmaMEff2 (0:ngllx-1))
    allocate(Enonli%specel(n)%SigmaEff1  (0:ngllx-1,6))
    allocate(Enonli%specel(n)%SigmaEff2  (0:ngllx-1,6))
    allocate(Enonli%specel(n)%Pore1      (0:ngllx-1))
    allocate(Enonli%specel(n)%Pore2      (0:ngllx-1))    
    allocate(Enonli%specel(n)%tauOCTO    (0:ngllx-1))    
    allocate(Enonli%specel(n)%gamOCTO    (0:ngllx-1))    
    allocate(Enonli%specel(n)%Sb         (0:ngllx-1))       !!!
    allocate(Enonli%specel(n)%aktifsur   (0:ngllx-1))       !!!



    Enonli%specel(n)%hexaSigma          = 0.0 
    Enonli%specel(n)%hexaS2             = 0.0
    Enonli%specel(n)%hexaS1             = 0.0
    Enonli%specel(n)%dhexaSigma         = 0.0
    Enonli%specel(n)%dhexaS             = 0.0
    Enonli%specel(n)%qStress1           = 0.0
    Enonli%specel(n)%qStress2           = 0.0
    Enonli%specel(n)%p0Stress           = 0.0
    Enonli%specel(n)%Gm0                = 0.0
    Enonli%specel(n)%S0                 = 1.0
    Enonli%specel(n)%S                  = 1.0
    Enonli%specel(n)%gammaref           = 0.0
    Enonli%specel(n)%eliflogic          = .True.
    Enonli%specel(n)%Wn                 = 0.0
    Enonli%specel(n)%r0                 = 0.0
    Enonli%specel(n)%r1                 = 0.0
    Enonli%specel(n)%S0old              = 1.0
    Enonli%specel(n)%Sold               = 1.0
    Enonli%specel(n)%Ws1                = 0.0
    Enonli%specel(n)%Ws2                = 0.0
    Enonli%specel(n)%SigmaP             = 0.0
    Enonli%specel(n)%Gmod1              = 0.0
    Enonli%specel(n)%Gmod2              = 0.0
    Enonli%specel(n)%Emod1              = 0.0
    Enonli%specel(n)%Emod2              = 0.0    
    Enonli%specel(n)%Kmod1              = 0.0
    Enonli%specel(n)%Kmod2              = 0.0  
    Enonli%specel(n)%gammaSAT           = 0.0
    Enonli%specel(n)%RSAT               = 0.0
    Enonli%specel(n)%CNinvSAT           = 0.0
    Enonli%specel(n)%SigmaMEff1         = 0.0
    Enonli%specel(n)%SigmaMEff2         = 0.0
    Enonli%specel(n)%SigmaEff1          = 0.0
    Enonli%specel(n)%SigmaEff2          = 0.0
    Enonli%specel(n)%Pore1              = 0.0
    Enonli%specel(n)%Pore2              = 0.0  
    Enonli%specel(n)%tauOCTO            = 0.0  
    Enonli%specel(n)%gamOCTO            = 0.0  
    Enonli%specel(n)%Sb                 = 0.0
    Enonli%specel(n)%aktifsur           = 0


    ! IWAN 1D3C
    allocate(Enonli%specel(n)%hexaE   (0:ngllx-1,6,6))
    allocate(Enonli%specel(n)%hexaEd  (0:ngllx-1,6,6))
    allocate(Enonli%specel(n)%plastSa1(0:ngllx-1,Enonli%Nspr,6))
    allocate(Enonli%specel(n)%plastSa2(0:ngllx-1,Enonli%Nspr,6))
    allocate(Enonli%specel(n)%plastF1 (0:ngllx-1,Enonli%Nspr))
    allocate(Enonli%specel(n)%plastF2 (0:ngllx-1,Enonli%Nspr))
    allocate(Enonli%specel(n)%dplastF1(0:ngllx-1,Enonli%Nspr))
    allocate(Enonli%specel(n)%dplastF2(0:ngllx-1,Enonli%Nspr))

    Enonli%specel(n)%hexaE    =   0.0
    Enonli%specel(n)%hexaEd   =   0.0
    Enonli%specel(n)%plastSa1 =   0.0
    Enonli%specel(n)%plastSa2 =   0.0
    Enonli%specel(n)%plastF1  =   0.0
    Enonli%specel(n)%plastF2  =   0.0
    Enonli%specel(n)%dplastF1 =   0.0
    Enonli%specel(n)%dplastF2 =   0.0 
enddo

return

end subroutine allocate_domain