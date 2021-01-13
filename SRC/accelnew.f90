subroutine accelnew(Tdomain,TimeP,Source,it,Evisco,Ewave,Enonli)

use sdomain 
use stimeparam
use ssource
use viscopara
use inciwave
use nonlinear

implicit none

type (domain),  intent (INOUT),  target ::  Tdomain
type (time),    intent (IN),     target ::  TimeP
type (sour),    intent (IN)             ::  Source
integer,        intent (IN)             ::  it
type(visco),    intent(INOUT),   target ::  Evisco
type(inci),     intent(IN),      target ::  Ewave
type(nonli),    intent(INOUT),   target ::  Enonli

integer         ::  ns
integer         ::  ncc 
integer         ::  i 
integer         ::  n 
integer         ::  itermax 
integer         ::  iter 
integer         ::  nd 
integer         ::  id 
integer         ::  nl1
integer         ::  nl2 
integer         ::  i_domain 
integer         ::  j
integer         ::  nw 
integer         ::  xx 
real            ::  Dt 
real            ::  Dt2
real            ::  alpha
real            ::  coef1
real            ::  coef2
real            ::  timealpha
real            ::  error
real            ::  errorinf
real            ::  vb
real            ::  F 
real            ::  F1
real            ::  F2
real            ::  rtime

integer, pointer          ::  ngllx 
integer, pointer          ::  nel 
integer, pointer          ::  nnode 
integer, pointer          ::  ngllxb  

real, dimension(:), allocatable     ::  s0
real, dimension(:), allocatable     ::  s1 
real, dimension(:), allocatable     ::  s2 
real, dimension(:), allocatable     ::  Vxloc
real, dimension(:), allocatable     ::  elVxloc
real, dimension(:), allocatable     ::  els1
real, dimension(:), allocatable     ::  els2
real, dimension(:), allocatable     ::  els0
real, dimension(:), allocatable     ::  zipto
real, dimension(:,:), allocatable   ::  HT 
real, dimension(:,:), allocatable   ::  H 

real, dimension(:), pointer         ::  A 
real, dimension(:), pointer         ::  elStress
real, dimension(:), pointer         ::  elMidDispl
real, dimension(:), pointer         ::  LocAcc
real, dimension(:), pointer         ::  DumpSx3

real, dimension(:,:), pointer       ::  LocCalc
real, dimension(:,:), pointer       ::  LocCalc0
real, dimension(:,:), pointer       ::  LocV0
real, dimension(:,:), pointer       ::  LocField
real, dimension(:,:), pointer       ::  LocFor
real, dimension(:,:), pointer       ::  LocMidDispl
real, dimension(:,:), pointer       ::  zipt 
real, dimension(:,:), pointer       ::  LocStress



! Predictor-MultiCorrector Scheme 
! Time staggered formulation Stress-Velocity in PML and contiguous elements
! Standard Newmark scheme in the medium.

nel     =>    Tdomain%nelem

alpha   =     TimeP%alpha
itermax =     TimeP%itermax

error     = -1.0
errorinf  = -1.0


do n = 0,nel-1
  ! PREDICTOR PHASE inside the bulk
  Dt    = Tdomain%specel(n)%dt
  Dt2   = Tdomain%specel(n)%dt**2

  coef1 = TimeP%beta/TimeP%gamm
  coef2 = 1.0/TimeP%gamm
             
  LocStress   => Tdomain%specel(n)%Stress 
  LocMidDispl => Tdomain%specel(n)%MidDispl     
  LocCalc     => Tdomain%specel(n)%Veloc 
  LocCalc0    => Tdomain%specel(n)%Calc0
  LocField    => Tdomain%specel(n)%Accel 
  LocV0       => Tdomain%specel(n)%V0

  ngllx       => Tdomain%specel(n)%ngllx

  LocV0       =  LocCalc                                       ! Vel n (if iter>0)
  LocCalc0    =  LocCalc                                       ! Vel 0,n+1
  LocField    =  (1.0-coef2)*LocField                          ! Accel 0,n+1

  ngllx       => Tdomain%specel(n)%ngllx 

  i_domain    = Tdomain%specel(n)%NumDomain-1
  HT          = Tdomain%Sub_domain(i_domain)%hTprime
  
  A           => Tdomain%specel(n)%Acoeff                           


    LocMidDispl(:,4) = LocStress(:,4) -Enonli%specel(n)%hexaSigma(:,4)
    LocMidDispl(:,5) = LocStress(:,5) -Enonli%specel(n)%hexaSigma(:,5)
    LocMidDispl(:,6) = LocStress(:,6) -Enonli%specel(n)%hexaSigma(:,6)
enddo


do n=0,Tdomain%nbdirichlet-1
! CORRECTION FOR DIRICHLET CONDITIONS
  vb    = -1.0
  rtime = float(it)*TimeP%dt 
  if (Tdomain%dirichlet(n)=="L") then
    call velocdiric(Tdomain%paramdiric(0,0),Tdomain%paramdiric(0,1) &
                                     ,Tdomain%paramdiric(0,2),rtime,vb)
    nd = 0
    id = 0
  else if (Tdomain%dirichlet(0)=="R") then
    call velocdiric(Tdomain%paramdiric(0,0),Tdomain%paramdiric(0,1) &
                                     ,Tdomain%paramdiric(0,2),rtime,vb)
    nd = nel-1
    id = Tdomain%specel(nd)%ngllx-1 
  else
    call velocdiric(Tdomain%paramdiric(1,0),Tdomain%paramdiric(1,1) &
                                     ,Tdomain%paramdiric(1,2),rtime,vb)
    nd = nel-1
    id = Tdomain%specel(nd)%ngllx-1
  endif 
  Dt  = Tdomain%specel(nd)%dt
  Dt2 = Tdomain%specel(nd)%dt**2
  do xx = 1,3
    Tdomain%specel(nd)%Veloc(id,xx) = vb
    Tdomain%specel(nd)%Accel(id,xx) = Tdomain%specel(nd)%Accel(id,xx) + (coef2/Dt)&
                                        *(vb-Tdomain%specel(nd)%V0(id,xx))
  enddo
enddo




if (Tdomain%borehole) then
! BOREHOLE CDT
  ngllx       => Tdomain%specel(nel-1)%ngllx 
  Dt    = Tdomain%specel(nel-1)%dt

  coef1 = TimeP%beta/TimeP%gamm
  coef2 = 1.0/TimeP%gamm

  Tdomain%specel(nel-1)%Veloc(ngllx-1,:) = Tdomain%Vborehole
  Tdomain%specel(nel-1)%Accel(ngllx-1,:) = Tdomain%specel(nel-1)%Accel(ngllx-1,:)+ (coef2/Dt)&
                              *(Tdomain%specel(nel-1)%Veloc(ngllx-1,:)-Tdomain%specel(nel-1)%V0(ngllx-1,:))
endif




do iter=0,itermax-1
  ! INTERNAL FORCES
  do n = 0,nel-1
    do xx = 1,3
      ngllx => Tdomain%specel(n)%ngllx
      allocate(s0(0:ngllx-1)) 
      allocate(H(0:ngllx-1,0:ngllx-1))
      s0 =  0.0
      H  =  0.0

      i_domain = Tdomain%specel(n)%NumDomain-1
      H        = Tdomain%Sub_domain(i_domain)%hprime

      A        => Tdomain%specel(n)%Acoeff
      s0       =  A* Tdomain%specel(n)%MidDispl(:,xx+3)
      Tdomain%specel(n)%Forces(:,xx) = MATMUL(H,s0)
      deallocate(s0,H)
    enddo
  enddo

  ! EXTERNAL FORCES 
  if ( Tdomain%source ) then
    do ns = 0,Source%ine-1
      ncc = Source%Elem(ns)%nr
      ngllx => Tdomain%specel(ncc)%ngllx
      timealpha = TimeP%rtime+alpha*Tdomain%specel(ncc)%dt  
      do xx=1,3

        if (Source%kind == 1) then
        ! Directional Impulse 
          do i = 0,ngllx-1
            Tdomain%specel(ncc)%Forces(i,xx) = Tdomain%specel(ncc)%Forces(i,xx) +   &
                                                    CompSource (Source,timealpha)*  &
                                                    Source%Elem(ns)%Impulse(i)
          enddo
        else if (Source%kind == 2 ) then
        ! Explosive Source
          do i = 0,ngllx-1
            Tdomain%specel(ncc)%Forces(i,xx) = Tdomain%specel(ncc)%Forces(i,xx) +   &
                            CompSource(Source,timealpha)* Source%Elem(ns)%Explosion(i)
          enddo    
        endif
      enddo
    enddo
  endif



  ! COMMUNICATING FORCES between the elements at each node
  do i=0,Tdomain%nnode-1
    if (Tdomain%node(i)%type == "I") then
      nl1 = Tdomain%node(i)%elemconnect(0)      
      nl2 = Tdomain%node(i)%elemconnect(1)
      ngllx => Tdomain%specel(nl1)%ngllx
      do xx=1,3
        F1 = Tdomain%specel(nl1)%Forces(ngllx-1,xx)
        F2 = Tdomain%specel(nl2)%Forces(0,xx)
        F = F1+F2
        Tdomain%specel(nl1)%Forces(ngllx-1,xx) = F
        Tdomain%specel(nl2)%Forces(0,xx) = F
      enddo
    endif
  enddo



  ! SOLUTION and CORRECTION PHASES
  do n = 0,nel-1

    Dt = Tdomain%specel(n)%dt
    Dt2 = Tdomain%specel(n)%dt**2
    ngllx => Tdomain%specel(n)%ngllx

    LocCalc     => Tdomain%specel(n)%Veloc 
    LocCalc0    => Tdomain%specel(n)%Calc0
    LocField    => Tdomain%specel(n)%Accel 
    LocV0       => Tdomain%specel(n)%V0
    LocStress   => Tdomain%specel(n)%Stress 
    LocFor      => Tdomain%specel(n)%Forces
    
    do xx = 1,3
      if (xx .NE. 3) then
          LocCalc(:,xx)   = Tdomain%specel(n)%DumpVx1 * LocV0(:,xx) &
                        + Tdomain%specel(n)%DumpVx2 * Dt * LocFor(:,xx) 

          LocField(:,xx)  = LocField(:,xx) + coef2/Dt*(LocCalc(:,xx)-LocCalc0(:,xx)) 
      else
          LocCalc(:,xx)   = Tdomain%specel(n)%DumpVx1E * LocV0(:,xx) &
                        + Tdomain%specel(n)%DumpVx2E * Dt * LocFor(:,xx) 

          LocField(:,xx)  = LocField(:,xx) + coef2/Dt*(LocCalc(:,xx)-LocCalc0(:,xx)) 
      endif

      if (iter < itermax-1) LocCalc0 = LocCalc 
    enddo
  enddo



  ! CORRECTION FOR DIRICHLET CONDITIONS
  do n=0,Tdomain%nbdirichlet-1
    rtime = float(it)*TimeP%dt 
    vb = -1.0
    if (Tdomain%dirichlet(n)=="L") then
      call velocdiric(Tdomain%paramdiric(0,0),Tdomain%paramdiric(0,1) &
                                     ,Tdomain%paramdiric(0,2),rtime,vb)
      nd = 0
      id = 0
    else if (Tdomain%dirichlet(0)=="R") then
      call velocdiric(Tdomain%paramdiric(0,0),Tdomain%paramdiric(0,1) &
                                     ,Tdomain%paramdiric(0,2),rtime,vb)
      nd = nel-1
      id = Tdomain%specel(nd)%ngllx-1 
    else
      call velocdiric(Tdomain%paramdiric(1,0),Tdomain%paramdiric(1,1) &
                                     ,Tdomain%paramdiric(1,2),rtime,vb)
      nd = nel-1
      id = Tdomain%specel(nd)%ngllx-1
    endif
    Dt  = Tdomain%specel(nd)%dt
    Dt2 = Tdomain%specel(nd)%dt**2

    do xx = 1,3    
      Tdomain%specel(nd)%Accel(id,xx) = Tdomain%specel(nd)%Accel(id,xx) + (coef2/Dt)&
                                          *(vb-Tdomain%specel(nd)%Veloc(id,xx))
      Tdomain%specel(nd)%Veloc(id,xx) = vb
      if (iter < itermax-1)  Tdomain%specel(nd)%Calc0(id,xx) = Tdomain%specel(nd)%Veloc(id,xx)
    enddo
  enddo


  if (Tdomain%borehole) then
  ! BOREHOLE CDT
    ngllx       => Tdomain%specel(nel-1)%ngllx 
    Dt    = Tdomain%specel(nel-1)%dt

    coef1 = TimeP%beta/TimeP%gamm
    coef2 = 1.0/TimeP%gamm


    Tdomain%specel(nel-1)%Accel(ngllx-1,:) = Tdomain%specel(nel-1)%Accel(ngllx-1,:)+ (coef2/Dt)&
                                *(Tdomain%Vborehole-Tdomain%specel(nel-1)%Veloc(ngllx-1,:))

    Tdomain%specel(nel-1)%Veloc(ngllx-1,:) = Tdomain%Vborehole   

    if (iter < itermax-1)  Tdomain%specel(nel-1)%Calc0(ngllx-1,:) = Tdomain%specel(nel-1)%Veloc(ngllx-1,:)
  endif

enddo ! iter

return

end subroutine accelnew