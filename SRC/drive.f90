program  drive_sem

! eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee !
! This file is the main file of the 1D-3C code that follows the procedure outlined in the 		  !
! manual Oral (2017) and the physics detailed in PhD dissertation of Oral (2016). 				      !
!
! eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee !


use sdomain
use stimeparam
use ssource
use sreceivers
use logical_input
use viscopara
use inciwave
use nonlinear

implicit none

type (domain),          target :: Tdomain
type (time),            target :: TimeP
type (sour),            target :: Source
type (receiver_stream), target :: Station
type (visco),           target :: Evisco
type (inci),            target :: Ewave
type (nonli),           target :: Enonli


character(len = 20)   ::  sourcefile
character(len = 20)   ::  stationsfile
character(len = 20)   ::  matfile
character(len = 20)   ::  meshfile

integer               ::  n
integer               ::  i 
integer               ::  ismod
integer               ::  j
integer               ::  i_domain
integer               ::  i_soil
integer               ::  idom 
integer               ::  it 
integer               ::  tfactor
integer               ::  nummat
integer               ::  k 
integer               ::  posi
integer               ::  eldim

logical               ::  irun

real                  ::  vs 
real                  ::  mu
real                  ::  Dt 
real                  ::  ellen
real                  ::  nt 
real                  ::  coef1
real                  ::  coef2
real                  ::  veloc1
real                  ::  veloc2
real                  ::  minuto0
real                  ::  minuto1

integer, pointer      ::  nel 
integer, pointer      ::  npt 
integer, pointer      ::  ngllx 
integer , pointer      ::  ndom
integer, pointer      ::  nspr

real, dimension(:), allocatable     ::  GlobCoord
real, dimension(:), allocatable     ::  Energy
real, dimension(:), allocatable     ::  Kin_Energy 
real, dimension(:), allocatable     ::  Pot_Energy 
real, dimension(:), allocatable     ::  R 
real, dimension(:), allocatable     ::  H
real, dimension(:), allocatable     ::  CNinv
real, dimension(:), allocatable     ::  s0 
real, dimension(:), allocatable     ::  s1 
real, dimension(:), allocatable     ::  s2 
real, dimension(:), allocatable     ::  Vxloc 
real, dimension(:), allocatable     ::  elVxloc
real, dimension(:), allocatable     ::  els1
real, dimension(:), allocatable     ::  els2 
real, dimension(:), allocatable     ::  els0
real, dimension(:), allocatable     ::  zipto

real, dimension(:,:), allocatable   ::  hT

real, dimension(:),   pointer       ::  DumpSx3
real, dimension(:,:), pointer       ::  zipt 



! ##############  Beginning of the program  #######################

call get_tiempo(minuto0)


! Creating the outfiles directory         ! Added FB 06/10/15

call execute_command_line ('mkdir -p outputfiles/')



! MODEL PROPERTIES

write (*,*)   "Read input file for model parameters"
open (11, file = "input.spec", status = "old")

read (11,*) Tdomain%Title_simulation
read (11,*) TimeP%duration
read (11,*) meshfile
read (11,*) matfile
read (11,*) Tdomain%source
read (11,*) sourcefile
read (11,*) stationsfile
read (11,*) Ewave%imposource
read (11,*) Ewave%imposig
read (11,*) Ewave%coeff
read (11,*) TimeP%skips
read (11,*) savetra
read (11,*) TimeP%output
read (11,*) irun
read (11,*) Tdomain%rheovisco
read (11,*) Tdomain%rheononli
read (11,*) Tdomain%Epressmod
read (11,*) Tdomain%Ewlevel
close(11)


TimeP%Newmark = 1
TimeP%itermax = 1
TimeP%beta  = 0.5
TimeP%gamm    = 1.0
TimeP%alpha   = 0.5



! MESH PROPERTIES 

write (*,*)   "Define mesh properties"
call define_mesh (Tdomain,meshfile,matfile)    
write (*,*)   "Compute Gauss-Lobatto-Legendre weights and zeroes"
call compute_GLL (Tdomain)


! SHAPE FUNCTIONS

npt => Tdomain%npoints                    
allocate(GlobCoord (0:npt-1))         
Globcoord = 0.0


! INITIAL CHARACTERISTIC SOIL MODULI

call soilmoduli (Tdomain)				


! VISCOELASTIC PROPERTIES

call viscoe (Evisco,Tdomain)


ndom      =>    Tdomain%nsubdomain
nel       =>    Tdomain%nelem
nspr      =>    Enonli%Nspr


! PRESSURE-DEPENDENT MODEL PPTS & IAI MODEL PARAMETERS

if (Tdomain%Epressmod) &
call readpress(Enonli)


! NONLINEAR PROPERTIES

if (Tdomain%rheononli) &
call NLnopress(Tdomain,Enonli,Evisco,i,TimeP)



do n = 0,nel-1
  ngllx     =>   Tdomain%specel(n)%ngllx
  allocate(Tdomain%specel(n)%InvGrad(0:ngllx-1))     
  allocate(Tdomain%specel(n)%Jacob  (0:ngllx-1))
enddo

write (*,*)   "Computing shape functions"
call shape2(TDomain,npt- 1,GlobCoord)
write (*,*)   "Computing derivatives of shape functions"
call deriv2(TDomain)


! Verify CFL

call TestCFL(Tdomain,TimeP,npt-1,GlobCoord)


! SOURCE PARAMETERS AND LOCATION

if (Tdomain%source) then
  write (*,*)   "Source parameters and location"
  call SourceParameters(Source,sourcefile)
  call SourcePosition(Tdomain,Source,GlobCoord,npt)
endif

TimeP%dt    =   Tdomain%dtmaximum
TimeP%nt    =   int(TimeP%Duration/ TimeP%dt)


! RECEIVER PARAMETERS AND LOCATION
eldim = int(TimeP%nt/ TimeP%output) +1 


if (savetra ) then
  write (*,*) "Receivers parameters and location"
  call ReceiverParam (Station,stationsfile)
  call ReceiverPosition(Tdomain,Station,GlobCoord,npt)
endif


call allocate_domain (Tdomain,Evisco,Enonli,Station,eldim)
 write(*,*) "Define Internal-External forces coefficients for Newmark scheme"
call define_Arrays (Tdomain)

if (saveenergy)   then
  allocate(Energy    (0:TimeP%nt-1))
  allocate(Kin_Energy(0:TimeP%nt-1))
  allocate(Pot_Energy(0:TimeP%nt-1))
else
  do n = 0,nel-1
    deallocate(Tdomain%specel(n)%Density)
  enddo
endif



! NEW 1D3C INITIAL CONDITIONS

call InitialCdts (Tdomain, Enonli, npt-1, GlobCoord, Evisco)


! IMPOSED SOURCE FILE

if (Ewave%imposource) &
call incidentw(Ewave,TimeP%dt)

 
if (irun) then
! ITERATION PROCESS NEWMARK SCHEME

print *,    "Entry Newmark scheme"
posi = 0

do i = 1,TimeP%nt
  TimeP%rtime = float(i-1)*TimeP%dt
     
  
  ! Attributing the input
  if (Ewave%imposource  .AND. Ewave%tfactor .ne. 1)   then
    if ( MOD((i-1),Ewave%tfactor) == 0 )    posi = posi+1
  else if ( Ewave%tfactor == 1) then
    posi = posi+ 1
  endif  
  
  ! Strain and/or stress calculation
  call sigeps(Tdomain,TimeP,Source,i,Evisco,Ewave,posi,Enonli)


  ! Stress for nonlinearity 1D-3C IWAN
  if (Tdomain%rheononli) &
  call Iwan3c(Tdomain,Enonli,i,TimeP,Evisco)
  

  ! IAI model
  if (Tdomain%Epressmod) then  
    call ShearWork(Tdomain,Enonli,i,TimeP,npt-1,Globcoord)
  endif  


  ! Final solution   
  call accelnew(Tdomain,TimeP,Source,i,Evisco,Ewave,Enonli)  

  ismod = mod(i,TimeP%skips)

  if (ismod == 0)   then
    print *, 'iter = ',i
    print *, 'elapsed time = ', i*TimeP%dt
  endif


  ! Storing traces
  if (savetra) call savetrace (Tdomain,Enonli,Station,i,TimeP)


enddo
close(58)
print *,"*************" 
print *,"Out of Newmark scheme"


! Unpack tracefile

if (savetra) then
  print *,"Saving traces"
  call dumptrace (Station,TimeP, Tdomain)
endif


endif ! irun


! Total elapsed time

call get_tiempo(minuto1)

write(*,*)
write(*,*) 'Elapsed Time (min) = ', minuto1-minuto0
write(*,*)

stop
end program 
