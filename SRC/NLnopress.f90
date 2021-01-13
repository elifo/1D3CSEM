subroutine NLnopress(Tdomain,Enonli,Evisco,temps,TimeP)

use sdomain
use nonlinear
use viscopara
use stimeparam


implicit none


type (domain),intent(INOUT),  target :: Tdomain
type (nonli), intent(INOUT),  target :: Enonli
type (visco), intent(IN),     target :: Evisco
integer,      intent(IN),     target :: temps 
type (time),  intent (IN),    target :: TimeP   


integer               ::  i
integer               ::  i_domain
integer               ::  idom 
integer               ::  i_soil
integer               ::  n
integer               ::  nummat
integer               ::  numsat
integer               ::  k

real                  ::  mu 
real                  ::  gammaref
real                  ::  x0
real                  ::  xu
real                  ::  dx
real                  ::  sumdummy
real                  ::  Gmod    

! integer, pointer      ::  ndom
integer, pointer      :: ndom

integer, pointer      ::  ngllx  
integer, pointer      ::  nel   
integer, pointer      ::  nspr  

real, dimension(:), allocatable     ::  R 
real, dimension(:), allocatable     ::  H
real, dimension(:), allocatable     ::  CNinv
real, dimension(:), pointer         ::  m1
real, dimension(:), allocatable     ::  strEA, bbcEA




ndom      =>    Tdomain%nsubdomain
nel       =>    Tdomain%nelem
nspr      =>    Enonli%Nspr

 

call get_G_info (Enonli)


allocate(Enonli%H(nspr,ndom))
allocate(Enonli%R(nspr,ndom))
allocate(Enonli%CNinv(nspr-1,ndom))

Enonli%H      =   0.0
Enonli%R      =   0.0
Enonli%CNinv  =   0.0

  


! Pressure-independent model (gref given )

if (.NOT. Tdomain%Epressmod)  then

  do i_domain = 0,ndom-1
 
    idom    =   i_domain+ 1
    i_soil  =   Tdomain%Sub_Domain(i_domain)%nonsol

    if ( .not. Tdomain%rheovisco ) then           ! Elastic Modulus
      mu = Tdomain%Sub_Domain(i_domain)%Mu
    else                                          ! Viscoelastic Modulus
      do n = 0,nel-1
        nummat = Tdomain%specel(n)%Numdomain- 1
        if (nummat == i_domain) then 
          mu = Evisco%specel(n)%visM(0)           ! Given identical Mu for each GLL point 
        exit
        endif
      enddo
      print*, "MU RATE (viscoelastic/elastic)=", mu/Tdomain%Sub_Domain(i_domain)%Mu 
    endif


    allocate(R(nspr))
    allocate(H(nspr))
    allocate(CNinv(nspr- 1))
    R      =  0.0
    H      =  0.0
    CNinv  =  0.0


    if ( Enonli%modeltype )   then
      open(33, file = 'outputfiles/SOILHYPER', status = 'unknown') 
      call nonline_curve_hyper(Tdomain,Enonli,mu,i_soil,R,H,CNinv)
    else
      open(33, file = 'outputfiles/SOILEXP', status = 'unknown') 
      call nonline_curve(Tdomain,Enonli,mu,i_soil,R,H,CNinv)
    endif

    Enonli%R     (:,idom)  = R
    Enonli%H     (:,idom)  = H
    Enonli%CNinv (:,idom)  = CNinv
    deallocate(R,H,CNinv)

    gammaref  =  Enonli%gammaref(i_soil)

  enddo
endif





end subroutine NLnopress
