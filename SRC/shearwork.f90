subroutine ShearWork(Tdomain,Enonli, it,TimeP,nm,GlobCoord)

use   sdomain
use   nonlinear 
use   stimeparam

implicit none

type(nonli),  intent(INOUT),  target    :: Enonli
type(domain), intent(INOUT),  target    :: Tdomain
integer,      intent(IN),     target    :: it
type (time),  intent (IN),    target    :: TimeP

integer,        intent (IN)             :: nm
real,           intent (IN)             :: GlobCoord(0:nm)


integer       ::  i 
integer       ::  n 
integer       ::  numsat
integer       ::  nummat
integer       ::  j 
integer       ::  nl1 
integer       ::  nl2

real          ::  yy1(6)
real          ::  depsv
real          ::  CoR
real          ::  dWse
real          ::  dWst
real          ::  dWs 
real          ::  w
real          ::  rr2 
real          ::  mm3
real          ::  rr3
real          ::  S2 
real          ::  Sparam
real          ::  sigmto
real          ::  WT

integer, pointer    ::  nel
integer, pointer    ::  ngllx 


real, dimension(:), pointer   ::  Cohesion
real, dimension(:), pointer   ::  poreold
real, dimension(:), pointer   ::  w1
real, dimension(:), pointer   ::  m1 
real, dimension(:), pointer   ::  m2
real, dimension(:), pointer   ::  c1 
real, dimension(:), pointer   ::  Wsold
real, dimension(:), pointer   ::  Sold
real, dimension(:), pointer   ::  S0old
real, dimension(:), pointer   ::  G0EAold
real, dimension(:), pointer   ::  G0EA
real, dimension(:), pointer   ::  Gm0
real, dimension(:), pointer   ::  Peff0
real, dimension(:), pointer   ::  sigmefold
real, dimension(:), pointer   ::  Told
real, dimension(:), pointer   ::  T 
real, dimension(:), pointer   ::  Q01
real, dimension(:), pointer   ::  Q02
real, dimension(:), pointer   ::  r 
real, dimension(:), pointer   ::  Ws 
real, dimension(:), pointer   ::  Wn
real, dimension(:), pointer   ::  S0
real, dimension(:), pointer   ::  p1
real, dimension(:), pointer   ::  p2
real, dimension(:), pointer   ::  ss1 
real, dimension(:), pointer   ::  S
real, dimension(:), pointer   ::  pore 
real, dimension(:), pointer   ::  sigmef
real, dimension(:), pointer   ::  tauOCTO
real, dimension(:), pointer   ::  gamOCTO


real, dimension(:,:), pointer     ::  sigma 
real, dimension(:,:), pointer     ::  eps 
real, dimension(:,:), pointer     ::  epsold
real, dimension(:,:), pointer     ::  sigmap  
real, dimension(:,:), pointer     ::  sigef 
real, dimension(:,:), pointer     ::  sigefold   	
real, dimension(:,:), pointer     ::  dplastic

real, dimension(:), pointer   ::  Sb
real          ::  ssl


integer :: nspr, k
real    :: sumdummy, gammaref
real, dimension(:,:), pointer   :: gammaSAT
real, dimension(:,:), pointer   :: RSAT
real, dimension(:,:), pointer   :: CNinvSAT



yy1(1)  = 1.0
yy1(2)  = 1.0
yy1(6)  = 1.0
yy1(3)  = 0.0
yy1(4)  = 0.0
yy1(5)  = 0.0


nel     =>    Tdomain%nelem

c1      =>    Enonli%SATc1 
m1      =>    Enonli%SATm1
m2      =>    Enonli%SATm2
w1      =>    Enonli%SATw1
p1      =>    Enonli%SATp1
p2      =>    Enonli%SATp2
ss1     =>    Enonli%SATs1 
Cohesion =>   Enonli%Cohesion


do n = 0,nel-1
  ngllx     =>  Tdomain%specel(n)%ngllx 
  nummat    =   Tdomain%specel(n)%Numdomain- 1 
  numsat    =   Tdomain%Sub_domain(nummat)%satsol

  sigma     =>  Tdomain%specel(n)%Stress
  eps       =>  Tdomain%specel(n)%gamma2
  epsold    =>  Tdomain%specel(n)%gamma1
  sigmap    =>  Enonli%specel(n)%SigmaP
  poreold   =>  Enonli%specel(n)%Pore1

  Wsold     =>  Enonli%specel(n)%Ws1
  Ws        =>  Enonli%specel(n)%Ws2
  Wn        =>  Enonli%specel(n)%Wn
  Sold      =>  Enonli%specel(n)%Sold
  S0old     =>  Enonli%specel(n)%S0old
  G0EAold   =>  Enonli%specel(n)%Gmod1
  G0EA      =>  Enonli%specel(n)%Gmod2
  Gm0       =>  Enonli%specel(n)%Gm0
  Peff0     =>  Enonli%specel(n)%p0Stress
  sigmefold =>  Enonli%specel(n)%SigmaMEff1
  Told      =>  Enonli%specel(n)%qStress1
  T         =>  Enonli%specel(n)%qStress2
  r         =>  Enonli%specel(n)%r1
  S0        =>  Enonli%specel(n)%S0
  S         =>  Enonli%specel(n)%S
  pore      =>  Enonli%specel(n)%Pore2
  sigefold  =>  Enonli%specel(n)%SigmaEff1
  sigef     =>  Enonli%specel(n)%SigmaEff2
  sigmef    =>  Enonli%specel(n)%SigmaMEff2
  tauOCTO   =>  Enonli%specel(n)%tauOCTO
  gamOCTO   =>  Enonli%specel(n)%gamOCTO


  dplastic  =>  Tdomain%specel(n)%dplastic

  Sb        =>  Enonli%specel(n)%Sb

  gammaSAT  =>  Enonli%specel(n)%gammaSAT
  RSAT      =>  Enonli%specel(n)%RSAT
  CNinvSAT  =>  Enonli%specel(n)%CNinvSAT  


	if ( Tdomain%Sub_domain(nummat)%Eeffective) then

    do i = 0,ngllx-1

      ! Principal stress
      call neoprincipals(sigma(i,:),sigmap(i,:))


      ! Deviatoric stress
      T(i)  = (sigmap(i,1)- sigmap(i,3))/ 2.0


      ! Plastic shear work increment
      dWs  = 0.0
      do j = 1,6
        dWs  = dWs+ (sigma(i,j)- yy1(j)*(sigma(i,1)+sigma(i,2)+sigma(i,6))/3.0)* dplastic(i,j)
      enddo



      ! Previous deviatoric ratio
      r(i)  = Told(i)/ Peff0(i)     

      ssl = 0.0
      ssl = 0.4+ (Sb(i)- 0.4)*S0old(i)/Sb(i)



      ! Correction of shear work
      if (Sold(i) .ge. ssl) then
        if (r(i)/S0old(i) .le. 0.67*m2(numsat)) then
          CoR = 1.0
        else
          CoR = (m1(numsat)-r(i)/Sold(i))/(m1(numsat)-0.67*m2(numsat))
        endif
      else
        if (r(i) .le. ssl*0.67*m2(numsat)) then
          CoR = 1.0
        else
          CoR = (ssl*m1(numsat) -r(i))/ (ssl*(m1(numsat)-0.67*m2(numsat))) 
        endif
      endif

      if (dWs  .gt. 0.0)&
      dWs = CoR* dWs


!!!
      dWs = max(0.0, dWs)




      Ws(i) = Wsold(i)+ dWs

      ! Normalized shear work
      w     = Ws(i)/ Wn(i)


      ! Variable S0
      if (w   .le.  0.0)    then
        S0(i)   =   1.0
      elseif  (w   .le.  w1(numsat))   then
        S0(i)   =   1.0- 0.6* ((w/w1(numsat))**p1(numsat))
      elseif  (w   .gt.  w1(numsat))   then
        S0(i)   =   (0.4- ss1(numsat))* ((w1(numsat)/w)**p2(numsat))+ ss1(numsat)
      endif


      ! Actual deviatoric ratio
      r(i)    =   T(i)/ Peff0(i)

      rr2     =   m2(numsat)* S0(i)
      mm3     =   m2(numsat)* 0.67
      rr3     =   mm3* S0(i)
      S2      =   S0(i)- (rr2-rr3)/ m1(numsat)

      ! VARIABLE S
      if (r(i) .le. rr3) then
        S(i)  = S0(i)
      else
        S(i)  = S2+ sqrt( (S0(i)- S2)**2.0+ ((r(i)- rr3)/ m1(numsat))**2.0 )
      endif


      !!! IF SOIL (GLL NODE) IS ABOVE WATER TABLE 
      WT = GlobCoord(Tdomain%specel(Tdomain%Ewlevel-1)%iglobnum(0))
      if (WT  >  GlobCoord(Tdomain%specel(n)%iglobnum(i))  ) then
        S(i)  = 1.0
        S0(i) = 1.0
      endif
      !

      ! Pore pressure & Effective stress
      pore(i)    = Peff0(i)* (1.0- S(i))
      do j = 1,6
        sigef(i,j) = sigma(i,j)- pore(i)*yy1(j)    
      enddo


      ! Updating saturation parameters
      poreold  = pore
      Told     = T
      Wsold    = Ws
      Sold     = S
      S0old    = S0 
      sigefold = sigef
      G0EAold  = G0EA


      ! reference strain
      gammaref = (m1(numsat)* Peff0(i) + cohesion(numsat)* Enonli%cos_phif(numsat))/ Gm0(i)  


      Enonli%specel(n)%gammaref(i) = gammaref

      nspr = Enonli%Nspr
  
      call neoupdate(S0old(i),Sb(i),Peff0(i),m1(numsat),&
            m2(numsat),Sold(i),gammaref,G0EA(i),nspr,gammaSAT(:,i),RSAT(:,i))

      sumdummy = 0.0
      do k = 1,nspr- 1
        CNinvSAT(k,i) = ((gammaSAT(k+ 1,i)/2.0- gammaSAT(k,i)/2.0)/(RSAT(k+ 1,i)- RSAT(k,i)))- 0.5/ G0EA(i)-sumdummy
        sumdummy      = sumdummy+ CNinvSAT(k,i)
      enddo

  enddo
  endif 
enddo 

end subroutine ShearWork
!
!
!

subroutine neoprincipals(sigma,sigmap)

implicit none

real,    intent(IN)      :: sigma(6)
real,    intent(INOUT)   :: sigmap(3)

integer  :: i
integer  :: j
integer  :: l

real     :: y(3) 
real     :: a
real     :: b
real     :: c
real     :: d
real     :: delta
real     :: gamma
real     :: TG
real     :: x1
real     :: x2
real     :: x3
real     :: x4
real     :: x5
real     :: x




if (sigma(3) == 0.0 .AND. sigma(4)== 0.0 .AND. sigma(5) == 0.0)   then
   y(1) = sigma(1)
   y(2) = sigma(2)
   y(3) = sigma(6)

do i = 1,2     
   l = i
   do j = i,3
      if (y(l) .gt. y(j) )   l = j
   enddo

   if(l.NE.i) then
      TG   = y(i)
      y(i) = y(l)
      y(l) = TG
   end if
end do 

do i = 1,3
   sigmap(i) = y(4-i)
end do

RETURN

endif



! Coefficients of the eqn of 3rd power
a = -1.0
b = sigma(1)+ sigma(2)+ sigma(6)
 c = sigma(3)**2.0 + sigma(4)**2.0+ sigma(5)**2.0 &
    - sigma(1)*sigma(2)&
    - sigma(2)*sigma(6)&
    - sigma(6)*sigma(1)
d = sigma(1)*sigma(2)*sigma(6) &
  + 2.0*sigma(3)*sigma(4)*sigma(5) &
  - sigma(1)*sigma(4)**2.0 &
  - sigma(2)*sigma(5)**2.0 &
  - sigma(6)*sigma(3)**2.0






delta = b**2.0- 3.0*a*c
if (delta .ne. 0.0) then
    gamma = (9.0*a*b*c- 2.0*b**3.0- 27.0*(a**2.0)*d)/(2.0* sqrt(abs(delta**3.0)))
endif

if (delta > 0.0) then
   if (abs(gamma) .le. 1.0) then
      x1 = (2.0* sqrt(delta)*(cos(acos(gamma)/3.0))-b)/(3.0*a)
      x2 = (2.0* sqrt(delta)* cos(acos(gamma)/3.0 - 2.* 4.0* atan(1.0)/3.0)-b)/(3.0*a)
      x3 = (2.0* sqrt(delta)* cos(acos(gamma)/3.0 + 2.* 4.0* atan(1.0)/3.0)-b)/(3.0*a)
   else
      RETURN
   endif
else if (delta .ge. 0.0) then
   RETURN
endif

y(1) = x1
y(2) = x2
y(3) = x3

do i = 1,2     
   l = i
   do j = i,3
      if (y(l) .gt. y(j) )   l = j
   enddo

   if(l.NE.i) then
      TG   = y(i)
      y(i) = y(l)
      y(l) = TG
   end if
end do 

do i = 1,3
   sigmap(i) = y(4-i)
end do

RETURN

end subroutine neoprincipals
!
!
!

subroutine neoupdate(S0old,Sb,Peff0,mm1,mm2,Sold,gammaref,G0EA,nspr,strEA,bbcEA)

implicit none

integer, intent(IN)        ::  nspr
real,    intent(IN)        ::  S0old
real,    intent(IN)        ::  Sb
real,    intent(IN)        ::  Peff0
real,    intent(IN)        ::  mm1 
real,    intent(IN)        ::  mm2
real,    intent(IN)        ::  Sold
real,    intent(IN)        ::  gammaref
real,    intent(INOUT)     ::  G0EA
real,    intent(INOUT)     ::  strEA(nspr)
real,    intent(INOUT)     ::  bbcEA(nspr)


integer     :: k
real        :: x0
real        :: xu
real        :: dx
real        :: delta
real        :: sigmax
real        :: maxeps
real        :: sumdummy , grefEA 



! [Iai et al. 1990 - Eqns 86-91]
delta   =   0.0
if (S0old .gt. Sb)   then
  sigmax    = Peff0*mm1*Sold
  grefEA    = gammaref
  G0EA      = sigmax/ grefEA         
else
  delta     = (mm1-mm2)*(Sb-S0old)*(0.4/Sb)* Peff0
  sigmax    = Peff0*mm1*Sold+ delta       
  grefEA    = gammaref/(S0old/ Sb)
  G0EA      = sigmax/ grefEA
endif

! Cn coefficients
x0      =   -6.0
xu      =   log10(0.10)     
dx      =   (xu-x0)/(nspr- 1)


do k = 1,nspr
  strEA(k)    =   10.0**(x0+ (k- 1)*dx)
  bbcEA(k)    =   G0EA* strEA(k)/ (1.0+ abs(strEA(k)/grefEA))
enddo

end subroutine neoupdate
