subroutine InitialCdts (Tdomain, Enonli, nm, GlobCoord,Evisco)

use sdomain
use nonlinear
use viscopara


implicit none


type(domain),   intent(IN),     target    :: Tdomain
type(nonli),    intent(INOUT),  target    :: Enonli
integer,        intent (IN)               :: nm
real,           intent (IN)               :: GlobCoord(0:nm)
type(visco) ,   intent (INOUT), target    :: Evisco



integer               ::  minz 
integer               ::  maxz
integer               ::  numsat
integer               ::  n
integer               ::  nummat
integer               ::  i
integer               ::  np
integer               ::  k
integer               ::  j
integer               ::  isoil
integer               ::  idom 
integer               ::  i_domain

real                  ::  rho 
real                  ::  a 
real                  ::  b 
real                  ::  c 
real                  ::  k0 
real                  ::  delta 
real                  ::  gammaref
real                  ::  thickness
real                  ::  x1
real                  ::  x2
real                  ::  mini
real                  ::  maxi
real                  ::  denomi
real                  ::  w0
real                  ::  Gmod
real                  ::  sumdummy
real                  ::  x0
real                  ::  xu
real                  ::  dx
real                  ::  WT

integer, pointer      ::  nel
integer, pointer      ::  ngllx  
integer , pointer      ::  ndom
integer      ::  nspr

real, dimension(:), allocatable   :: R 
real, dimension(:), allocatable   :: H
real, dimension(:), allocatable   :: CNinv
real, dimension(:), allocatable   :: strEA
real, dimension(:), allocatable   :: bbcEA

real, dimension(:), pointer       :: cohesion
real, dimension(:), pointer       :: m1 
real, dimension(:), pointer       :: m2
real, dimension(:), pointer       :: p1  
real, dimension(:), pointer       :: p2
real, dimension(:), pointer       :: ss1
real, dimension(:), pointer       :: w1
real, dimension(:), pointer       :: r0
real, dimension(:), pointer       :: S0old
real, dimension(:), pointer       :: qStress1
real, dimension(:), pointer       :: qStress2
real, dimension(:), pointer       :: Wn
real, dimension(:), pointer       :: T
real, dimension(:), pointer       :: Told  
real, dimension(:), pointer       :: Peff0
real, dimension(:), pointer       :: Gm0
real, dimension(:), pointer       :: S0
real, dimension(:), pointer       :: S 

real, dimension(:,:), pointer     :: sigmaini
real, dimension(:,:), pointer     :: S2
real, dimension(:,:), pointer     :: S1
real, dimension(:,:), pointer     :: dhexaSigma
real, dimension(:,:), pointer     :: dhexaS


real  :: m3
real  :: m4
real, dimension(:), pointer       :: Sb






nel     =>    Tdomain%nelem
! Initial stress computation

do n = 0,nel-1
  ngllx       =>    Tdomain%specel(n)%ngllx
  ndom        =>    Tdomain%nsubdomain

  nummat      =     Tdomain%specel(n)%Numdomain- 1
  rho         =     Tdomain%Sub_domain(nummat)%rho    

  sigmaini    =>    Enonli%specel(n)%hexaSigma


  ! This is gonna change if Dilatancy + Pressure-independent model

  if (n .ge. Tdomain%Ewlevel-1) &
  rho = Tdomain%Sub_domain(nummat)%rho- 1000.0

  if (.NOT. Tdomain%Epressmod) &
  rho = Tdomain%Sub_domain(nummat)%rho


  if (Tdomain%Epressmod) then
    numsat = Tdomain%Sub_domain(nummat)%satsol
    k0     = Enonli%Isoconso(numsat)
    print*, 'Element', n+1, 'K0', k0
  else
    k0     = 1.0
  endif 

  
  do i = 0,ngllx-1
    np =  Tdomain%specel(n)%iglobnum(i)
    if (i == 0 .AND. n == 0) then
      sigmaini(i,:) = 0.0
    else if (i == 0) then
      sigmaini(i,:) = Enonli%specel(n-1)%hexaSigma(ngllx-1,:)
    else
      sigmaini(i,1) = sigmaini(i-1,1)+ rho* k0* 9.8* (GlobCoord(np)-GlobCoord(np-1))
      sigmaini(i,2) = sigmaini(i-1,2)+ rho* k0* 9.8* (GlobCoord(np)-GlobCoord(np-1))
      sigmaini(i,6) = sigmaini(i-1,6)+ rho* 9.8* (GlobCoord(np)-GlobCoord(np-1))
    endif
  enddo 
enddo

! Nonzero surface node
Enonli%specel(0)%hexaSigma(0,:) = Enonli%specel(0)%hexaSigma(1,:)


! Midstress computation
Tdomain%Sub_domain(0)%midStress1   = 0.0
Tdomain%Sub_domain(0)%midStress2   = 0.0

do n = 0,ndom-1
  minz  = 0
  maxz  = 0
  rho   = Tdomain%Sub_domain(n)%rho

  j = 0
  do i = 0,nel-1
    nummat = Tdomain%specel(i)%Numdomain-1
    if (nummat .eq. n) then

      if (j .eq. 0)     minz = i
      if ( minz .gt. i) minz = i         
      if ( maxz .lt. i) maxz = i
      j = j+ 1
   
      if (i .ge. Tdomain%Ewlevel-1) &
      rho = Tdomain%Sub_domain(nummat)%rho- 1000.0

      if (.NOT. Tdomain%Epressmod) &
      rho = Tdomain%Sub_domain(nummat)%rho

    endif
  enddo

  minz      = Tdomain%specel(minz)%iglobnum(0)
  maxz      = Tdomain%specel(maxz)%iglobnum(ngllx- 1)
  thickness = GlobCoord(maxz)- GlobCoord(minz)

  Tdomain%Sub_domain(n)%midStress1 = rho* 9.8* thickness* (1.0+ 2.0* k0)/ 3.0

  ! IF WATER TABLE INSIDE THE DOMAIN
  WT = GlobCoord(Tdomain%specel(Tdomain%Ewlevel-1)%iglobnum(0))
  if (WT > GlobCoord(minz)  .AND.  WT < GlobCoord(maxz)) then
     Tdomain%Sub_domain(n)%midStress1 = Tdomain%Sub_domain(n)%midStress1+ &
                                    9.8* 1000.0* (WT-GlobCoord(minz))* (1.0+ 2.0* k0)/ 3.0          
  endif

  if (n .ne. 0)     &
  Tdomain%Sub_domain(n)%midStress1 =   Tdomain%Sub_domain(n)%midStress1 +&
                                          Tdomain%Sub_domain(n- 1)%midStress1
enddo                                              

do n = 0,ndom-1
  if (n == 0) then
    Tdomain%Sub_domain(n)%midStress2 =  Tdomain%Sub_domain(n)%midStress1/ 2.0   
  else
    Tdomain%Sub_domain(n)%midStress2 =  (Tdomain%Sub_domain(n)%midStress1+ &
                                          Tdomain%Sub_domain(n-1)%midStress1)/ 2.0  
  endif
enddo




! IF TOTAL STRESS ANALYSIS

if (.NOT. Tdomain%Epressmod) then  
  do n = 0,nel- 1
    ngllx       =>   Tdomain%specel(n)%ngllx

    sigmaini    =>   Enonli%specel(n)%hexaSigma
    S2          =>   Enonli%specel(n)%hexaS2
    S1          =>   Enonli%specel(n)%hexaS1
    T           =>   Enonli%specel(n)%qStress2
    Peff0       =>   Enonli%specel(n)%p0Stress
    Gm0         =>   Enonli%specel(n)%Gm0

    do i = 0,ngllx- 1
      T(i)      =   (sigmaini(i,6)- sigmaini(i,2))/2.0
      Peff0(i)  =   (sigmaini(i,1)+ sigmaini(i,2)+ sigmaini(i,6))/3.0
      Gm0(i)    =   Tdomain%specel(n)%Gmodulus

      ! Viscosity in
      if (Tdomain%rheovisco) then
        Gm0(i)  = Evisco%specel(n)%visM(i)
      endif

      S1(i,1)   =   sigmaini(i,1)- Peff0(i)
      S1(i,2)   =   sigmaini(i,2)- Peff0(i)
      S1(i,6)   =   sigmaini(i,6)- Peff0(i)
        
      do k = 1,6
        S1(i,k) =   max(0.001, S1(i,k))
      enddo
    enddo 

    if (Tdomain%specel(n)%PML)  sigmaini = 0.0
    
    Tdomain%specel(n)%Stress = sigmaini 
  enddo
endif


! IF PRESSURE-DEPENDENT MODEL
if (Tdomain%Epressmod) then
  do n = 0,nel-1
    ngllx       =>    Tdomain%specel(n)%ngllx
    nummat      =     Tdomain%specel(n)%Numdomain- 1
    numsat      =     Tdomain%Sub_domain(nummat)%satsol
    isoil       =     Tdomain%Sub_Domain(nummat)%nonsol

    sigmaini    =>    Enonli%specel(n)%hexaSigma
    Peff0       =>    Enonli%specel(n)%p0Stress
    Told        =>    Enonli%specel(n)%qStress1
    Gm0         =>    Enonli%specel(n)%Gm0
    Wn          =>    Enonli%specel(n)%Wn
    S1          =>    Enonli%specel(n)%hexaS1    

    r0          =>    Enonli%specel(n)%r0
    m1          =>    Enonli%SATm1
    m2          =>    Enonli%SATm2
    S0old       =>    Enonli%specel(n)%S0old
    ss1         =>    Enonli%SATs1
    w1          =>    Enonli%SATw1
    p1          =>    Enonli%SATp1
    p2          =>    Enonli%SATp2
    Sb          =>    Enonli%specel(n)%Sb

    cohesion    =>    Enonli%Cohesion


    do i = 0, ngllx-1
      Told  (i)   =   abs(sigmaini(i,6)- sigmaini(i,2))/ 2.0
      Peff0 (i)   =   (sigmaini(i,1)+ sigmaini(i,2)+ sigmaini(i,6))/ 3.0
      Gm0   (i)   =   Tdomain%specel(n)%Gmodulus &
                          * ((abs(Peff0(i)/Tdomain%Sub_domain(nummat)%midStress2))**0.5)

      ! Viscosity in
      if (Tdomain%rheovisco) then
         Gm0(i)  = Evisco%specel(n)%visM(i) &
                          * ((abs(Peff0(i)/Tdomain%Sub_domain(nummat)%midStress2))**0.5)
      endif



      S1(i,1)   =   sigmaini(i,1)- Peff0(i)
      S1(i,2)   =   sigmaini(i,2)- Peff0(i)
      S1(i,6)   =   sigmaini(i,6)- Peff0(i)
        
      do k = 1,6
        S1(i,k) =   max(0.001, S1(i,k))
      enddo

      Wn(i)   =   ((Peff0(i)* m1(numsat))**2.0) / (2.0* Gm0(i))

      ! Initial shear work
      r0(i)   =   Told(i)/ Peff0(i)


      ! After Fabian's NOAH code
      m4 = 0.0
      m3 = 0.67* m2(numsat)
      if (r0(i) .le. m3   .OR.  m1(numsat) == 0.0 ) then
        Enonli%specel(n)%Ws1(i) = 0.0
        S0old(i) = 1.0
        Sb(i)    = 0.4

      else
        m4 = 1.0- (m2(numsat)-m3)/m1(numsat)
        a  = m4*m4*m1(numsat)*m1(numsat)- m2(numsat)*m2(numsat)- 2.0*m3*m3+ 2.0*m2(numsat)*m3
        b  = 2.0*r0(i)*m3- 2.0*m1(numsat)*m1(numsat)*m4
        c  = m1(numsat)*m1(numsat)- r0(i)*r0(i)

       
        S0old(i) = (-b- sqrt(b*b- 4.0*a*c))/(2.0*a)

        if (S0old(i) .le. ss1(numsat)) &
        S0old(i) = ss1(numsat)

        if (S0old(i) .ge. 0.4) then
          Sb(i) = 0.4
        else
          Sb(i) = S0old(i)
        endif

        if (S0old(i) .ge. 0.4) then
          w0 = w1(numsat)* ((1.0-S0old(i))/0.6)**(1.0/p1(numsat))
        else
          if (S0old(i) .le. ss1(numsat)) then
            w0 = 1.0
          else
            w0 = w1(numsat)* ((S0old(i)-ss1(numsat))/(0.4-ss1(numsat)))**(-1.0/p2(numsat))
          endif         
        endif

        Enonli%specel(n)%Ws1(i) = w0* Wn(i)      

      endif



      if (.not.  Tdomain%Sub_domain(nummat)%Eeffective) then

        idom     = nummat+ 1

        gammaref = (m1(numsat)* Peff0(i) + cohesion(numsat)* Enonli%cos_phif(numsat))/ Gm0(i)   
		
        Enonli%specel(n)%gammaref(i) = gammaref
        nspr = Enonli%Nspr

        allocate(strEA(nspr))
        allocate(bbcEA(nspr))
        allocate(CNinv(nspr- 1))
        strEA = 0.0
        bbcEA = 0.0
        CNinv = 0.0

        ! Cn coefficients
        x0    = -6.0
        xu    = log10(0.10)     
        dx    = (xu-x0)/(nspr- 1)

        do k = 1,nspr
          strEA(k)    =   10.0**(x0+ (k- 1)*dx)   
          bbcEA(k)    =   Gm0(i)* strEA(k)/ (1.0+ abs(strEA(k)/gammaref))

          if (gammaref .lt. 1.e-15) &
          print*, "CHECK GAMMA_REFERENCE: NaN Error there"
        enddo

        sumdummy    = 0.0
        do k = 1,nspr- 1
          CNinv(k)  = ((strEA(k+ 1)/2.0- strEA(k)/2.0)/(bbcEA(k+ 1)- bbcEA(k)))- 0.5/ Gm0(i)-sumdummy
          sumdummy  = sumdummy+ CNinv(k)
        enddo

        Enonli%specel(n)%RSAT(:,i) = bbcEA
        Enonli%specel(n)%CNinvSAT(:,i) = CNinv

        deallocate(strEA,bbcEA,CNinv)

      endif

    enddo   
    write(99,*)

    !!!
    Enonli%specel(n)%Gmod2        =   Gm0
    Enonli%specel(n)%Gmod1        =   Gm0
    Enonli%specel(n)%SigmaMEff1   =   Peff0
    Enonli%specel(n)%SigmaEff1    =   sigmaini

    if (Tdomain%specel(n)%PML)  sigmaini = 0.0 
    Tdomain%specel(n)%Stress      =   sigmaini  

  enddo
endif

end subroutine InitialCdts
