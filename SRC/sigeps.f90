subroutine sigeps(Tdomain,TimeP,Source,it,Evisco,Ewave,posi,Enonli)

use sdomain 
use stimeparam
use ssource
use viscopara
use inciwave
use nonlinear


implicit none

type (domain),  intent (INOUT), target    ::  Tdomain
type (time),    intent (IN),    target    ::  TimeP
type (sour),    intent (IN)               ::  Source
integer, intent (IN)                      ::  it
type(visco), intent(INOUT),     target    ::  Evisco
type(inci), intent(INOUT),      target    ::  Ewave
integer, intent (IN)                      ::  posi
type(nonli),    intent(INOUT),  target    ::  Enonli


integer                   ::  n
integer                   ::  i_domain
integer                   ::  j 
integer                   ::  xx
integer                   ::  tfactor
integer                   ::  k 
integer                   ::  nummat
integer                   ::  interpos

real                      ::  Dt 
real                      ::  coef1 
real                      ::  coef2
real                      ::  interpo
real                      ::  yenihiz
real                      ::  veloc1
real                      ::  veloc2

integer, pointer                  ::  nel 
integer, pointer                  ::  npt 
integer, pointer                  ::  ngllx
integer, pointer                  ::  ndom 

real, dimension(:), allocatable   ::  ni 
real, dimension(:), allocatable   ::  lambda
real, dimension(:), allocatable   ::  mu 

real, dimension(:), allocatable   ::  s1
real, dimension(:), allocatable   ::  Vxloc 
real, dimension(:), allocatable   ::  elVxloc 
real, dimension(:), allocatable   ::  els1
real, dimension(:), allocatable   ::  els2
real, dimension(:), allocatable   ::  s2 
real, dimension(:), allocatable   ::  zipto
real, dimension(:,:), allocatable ::  HT 

real, dimension(:), pointer       ::  DumpSx3 => null()
real, dimension(:), pointer       ::  DumpSx3YZ=> null()
real, dimension(:), pointer       ::  DumpSx3ZZ=> null()
real, dimension(:), pointer       ::  DumpSx3VOL=> null()


real, dimension(:,:), pointer     ::  zipt=> null()
real, dimension(:,:), pointer     ::  ziptYZ=> null()
real, dimension(:,:), pointer     ::  ziptZZ=> null()
real, dimension(:,:), pointer     ::  ziptXX=> null()


real, dimension(:,:), pointer     ::  LocStress 
real, dimension(:,:), pointer     ::  Veloc 
real, dimension(:,:), pointer     ::  Accel

real, dimension(:), allocatable   ::  sig1
real, dimension(:), allocatable   ::  sig2
real  :: coeff



nel   =>   Tdomain%nelem

! STRAIN and/or STRESS CALCULATION

do n = 0,nel-1  
    
  ngllx       =>  Tdomain%specel(n)%ngllx
  LocStress   =>  Tdomain%specel(n)%Stress 
  Veloc       =>  Tdomain%specel(n)%Veloc
  Accel       =>  Tdomain%specel(n)%Accel 

  Dt      =   Tdomain%specel(n)%dt
  coef1   =   TimeP%beta/TimeP%gamm
  coef2   =   1.0/TimeP%gamm

  ! Updating previous strain
	Tdomain%specel(n)%gamma1 =  Tdomain%specel(n)%gamma2
	Tdomain%specel(n)%gamma2 =  0.0



  
  do xx = 1,3

    allocate(HT(0:ngllx-1,0:ngllx-1))
    allocate(s1(0:ngllx-1))
    allocate(s2(0:ngllx-1)) 
    allocate(Vxloc(0:ngllx-1))
    allocate(elVxloc(0:ngllx-1))
    allocate(els1(0:ngllx-1))
    allocate(els2(0:ngllx-1))  

    allocate(ni(0:ngllx-1))
    allocate(lambda(0:ngllx-1))
    allocate(mu(0:ngllx-1))

    allocate(sig1(0:ngllx-1))
    allocate(sig2(0:ngllx-1))




    ni = Tdomain%specel(n)%Ni 
    mu = Tdomain%specel(n)%Mu 
    lambda = -2.0* Tdomain%specel(n)%Mu+ Tdomain%specel(n)%lambda2mu


    s1      = 0.0
    s2      = 0.0 
    HT      = 0.0 
    Vxloc   = 0.0 
    elVxloc = 0.0
    els1    = 0.0 
    els2    = 0.0
 
    i_domain  =   Tdomain%specel(n)%NumDomain- 1
    HT        =   Tdomain%Sub_domain(i_domain)%hTprime               
    Vxloc     =   Veloc(:,xx)+ Dt*(0.5- coef1)*(1.0- coef2)* Accel(:,xx)



    ! ABSORBING LAYER (PML) USED WITH INCIDENT WAVE (of 1 element of PML)
    ! input velocity multiplied by 2 (shared by upside and downside)
    if (Tdomain%specel(nel-1)%PML) then
    Ewave%sposi = 1
        if (Ewave%imposource  .AND. Ewave%tfactor .ne. 1) then 
        ! Input signal implementation & interpolation                     
        elVxloc    =   0.0
        interpos   =   MOD((it- 1),Ewave%tfactor)
        interpo    =   float(interpos)/ float(Ewave%tfactor) 

        if (n == nel- Ewave%sposi-1  .AND. TimeP%rtime .le. Ewave%timeup)   then     

          veloc1           =   Ewave%VB(posi,xx)*Ewave%coeff
          veloc2           =   Ewave%VB(posi+1,xx)*Ewave%coeff               
          yenihiz          =   veloc1+ (interpo)*(veloc2- veloc1)
          elVxloc(ngllx-1)       =   yenihiz* 2.0    

          els1             =   MATMUL(hT,elVxloc)
          els2             =   Tdomain%specel(n)%InvGrad* els1
        else
          els2 = 0.0
        endif 
      else if (Ewave%tfactor == 1)    then
        elVxloc  =  0.0
        if ( n == nel- Ewave%sposi-1  .AND. TimeP%rtime .le. Ewave%timeup)   then

          yenihiz          =   Ewave%VB(posi,xx)* Ewave%coeff 
          elVxloc(ngllx-1) =   yenihiz* 2.0

          els1             =   MATMUL(hT,elVxloc)
          els2             =   Tdomain%specel(n)%InvGrad* els1
        else 
          els2 = 0.0
        endif 
      else if (.not. Ewave%imposource)    then
        els2 = 0.0
      endif  
    endif


    ! RIGID BOUNDARY 
    ! (NORMALLY INCIDENT WAVE IS DEFINED ON THE LAST POINT, SO NO MULTIPLICATION BY 2)
    ! For some other tests, it is useful to define it on different elements
    ! thus, Ewave%sposi is still kept here !

    if ( .NOT. Tdomain%specel(nel-1)%PML  .AND. .not. Tdomain%borehole) then
        Ewave%sposi = 0
        if (Ewave%imposource  .AND. Ewave%tfactor .ne. 1) then 
          ! Input signal implementation & interpolation                     
          elVxloc    =   0.0
          interpos   =   MOD((it- 1),Ewave%tfactor)
          interpo    =   float(interpos)/ float(Ewave%tfactor) 

         if (n == nel- 1- Ewave%sposi  .AND. TimeP%rtime .le. Ewave%timeup)   then        

            veloc1           =   Ewave%VB(posi,xx)*Ewave%coeff
            veloc2           =   Ewave%VB(posi+1,xx)*Ewave%coeff               
            yenihiz          =   veloc1+ (interpo)*(veloc2- veloc1)
            elVxloc(ngllx-1) =   yenihiz


            els1             =   MATMUL(hT,elVxloc)
            els2             =   Tdomain%specel(n)%InvGrad* els1
          else
            els2 = 0.0
          endif 
        else if (Ewave%tfactor == 1)    then
          elVxloc  =  0.0
          if ( n == nel- 1- Ewave%sposi  .AND. TimeP%rtime .le. Ewave%timeup)   then     

            yenihiz          =   Ewave%VB(posi,xx)* Ewave%coeff 
            elVxloc(ngllx-1) =   yenihiz

            els1             =   MATMUL(hT,elVxloc)
            els2             =   Tdomain%specel(n)%InvGrad* els1
          else 
            els2 = 0.0
          endif 
        else if (.not. Ewave%imposource)    then
          els2 = 0.0
        endif  
    endif


    if (Tdomain%borehole) then
    ! Borehole uses always the last point of the model (bottom)
    ! Ewave%sposi is neglected
        yenihiz = 0.0
        if (Ewave%imposource  .AND. Ewave%tfactor .ne. 1) then 
          ! Input signal implementation & interpolation                     
          interpos   =   MOD((it- 1),Ewave%tfactor)
          interpo    =   float(interpos)/ float(Ewave%tfactor) 
          if (n == nel- 1  .AND. TimeP%rtime .le. Ewave%timeup)   then        
            veloc1           =   Ewave%VB(posi,xx)*Ewave%coeff
            veloc2           =   Ewave%VB(posi+1,xx)*Ewave%coeff               
            yenihiz          =   veloc1+ (interpo)*(veloc2- veloc1)
          endif
        else if (Ewave%imposource .AND. Ewave%tfactor == 1)    then
          if ( n == nel- 1  .AND. TimeP%rtime .le. Ewave%timeup)   then     
            yenihiz          =   Ewave%VB(posi,xx)* Ewave%coeff 
          endif 
        endif 
        Tdomain%Vborehole(xx) = yenihiz 
        els2 = 0.0
    endif


    ! Internal strain rate 
    s1 = MATMUL(hT,Vxloc)
    s2 = Tdomain%specel(n)%InvGrad* s1

    ! Strain rate contribution by input signal
    s2 = s2 + els2


    ! PML CONDITION
    if ( Tdomain%specel(n)%PML ) then
      if (xx .NE. 3) then
          LocStress(:,xx+ 3) = Tdomain%specel(n)%DumpSx1* LocStress(:,xx+ 3)+ &           
                                Tdomain%specel(n)%DumpSx2* Dt* s2* Tdomain%specel(n)%Mu
      else
          LocStress(:,xx+3) = Tdomain%specel(n)%DumpSx1E* LocStress(:,xx+3)+ &           
                       Tdomain%specel(n)%DumpSx2E* Dt* s2* Tdomain%specel(n)%lambda2mu
      endif
    endif


    ! HERE IT BEGINS 
    if ( .not. Tdomain%specel(n)%PML ) then
      if (.not. Tdomain%rheovisco  .AND.  .not. Tdomain%rheononli) then

        ! ELASTICITY
        Tdomain%specel(n)%gamma2(:,xx+ 3) = s2

        if (xx .NE. 3) then
          LocStress(:,xx+ 3) = Tdomain%specel(n)%DumpSx1* LocStress(:,xx+ 3)+ &           
                                Tdomain%specel(n)%DumpSx2* Dt* s2* Tdomain%specel(n)%Mu
        else
          LocStress(:,1) = Tdomain%specel(n)%DumpSx1* LocStress(:,1)+ &           
                       Tdomain%specel(n)%DumpSx2* Dt* s2* lambda

          LocStress(:,2) = Tdomain%specel(n)%DumpSx1* LocStress(:,2)+ &           
                       Tdomain%specel(n)%DumpSx2* Dt* s2* lambda

          LocStress(:,6) = Tdomain%specel(n)%DumpSx1* LocStress(:,6)+ &           
                       Tdomain%specel(n)%DumpSx2* Dt* s2* Tdomain%specel(n)%lambda2mu
        endif

      else if ( Tdomain%rheovisco ) then
        ! Visco pointers & arrays
        DumpSx3   => Evisco%specel(n)%DumpSx3
        DumpSx3YZ => Evisco%specel(n)%DumpSx3YZ
        DumpSx3ZZ => Evisco%specel(n)%DumpSx3ZZ
        DumpSx3VOL=> Evisco%specel(n)%DumpSx3VOL

        zipt    => Evisco%specel(n)%viszipt
        ziptYZ  => Evisco%specel(n)%visziptYZ
        ziptZZ  => Evisco%specel(n)%visziptZZ
        ziptXX  => Evisco%specel(n)%visziptVOL

        allocate(zipto(0:ngllx-1)) 
        zipto   = 0.0  

        if (xx == 1) then
          zipto   = 0.0  
          ! Viscoelastic strain rate XZ
          do k = 0,ngllx-1
            do j = 1,8
              zipt(k,j) = zipt(k,j)* exp(-Dt/ Evisco%vistau(j))+ Evisco%specel(n)%visw(j)&
                            *(1.0- exp(-Dt/ Evisco%vistau(j)))*s2(k)
              zipto(k)  = zipto(k)+ zipt(k,j)
            enddo
            DumpSx3(k)  =  (s2(k)- zipto(k))
          enddo  
        else if (xx == 2) then
          zipto   = 0.0   
          ! Viscoelastic strain rate YZ
          do k = 0,ngllx-1
            do j = 1,8
              ziptYZ(k,j) = ziptYZ(k,j)* exp(-Dt/ Evisco%vistau(j))+ Evisco%specel(n)%visw(j)&
                              *(1.0- exp(-Dt/ Evisco%vistau(j)))*s2(k)
              zipto(k)  = zipto(k)+ ziptYZ(k,j)
            enddo
            DumpSx3YZ(k)  =  (s2(k)- zipto(k))  
          enddo 
        else

          zipto   = 0.0            

          do k = 0,ngllx-1

          !xx
          call MAT_VISLA_strain2(Evisco%specel(n)%visziptVOL(k,:), Dt, Evisco%specel(n)%viswp, Evisco%specel(n)%visw, &
                                0.0, s2(k), Evisco%specel(n)%visMp(k), Evisco%specel(n)%visM(k),sig1(k))
          !zz
          call MAT_VISLA_strain2(Evisco%specel(n)%visziptZZ(k,:), Dt, Evisco%specel(n)%viswp, Evisco%specel(n)%visw, &
                                s2(k), 0.0, Evisco%specel(n)%visMp(k), Evisco%specel(n)%visM(k),sig2(k))


          ! Used for visco-elastoplasticity
          coeff = 1.0 / (Evisco%specel(n)%visMp(k) *Evisco%specel(n)%visMp(k) &
                         - (Evisco%specel(n)%visMp(k) - 2d0 * Evisco%specel(n)%visM(k)) &
                         * (Evisco%specel(n)%visMp(k) - 2d0 * Evisco%specel(n)%visM(k)))

          ! Viscoelastic strain rate ZZ
          DumpSx3ZZ(k) = coeff * (Evisco%specel(n)%visMp(k) * sig2(k) & 
                               - (Evisco%specel(n)%visMp(k) - 2d0 * Evisco%specel(n)%visM(k)) * sig1(k))
          enddo

        endif

        
              

        if (.not. Tdomain%rheononli ) then
          ! VISCOELASTICITY
          Tdomain%specel(n)%gamma2(:,xx+ 3)    = s2 

          if (xx == 1) then !xz
            LocStress(:,xx+ 3) = Tdomain%specel(n)%DumpSx1* LocStress(:,xx+ 3) + &
                                    Tdomain%specel(n)%DumpSx2* Dt* Evisco%specel(n)%visM* DumpSx3 

          else if (xx == 2) then !yz
            LocStress(:,xx+ 3) = Tdomain%specel(n)%DumpSx1* LocStress(:,xx+ 3) + &
                                    Tdomain%specel(n)%DumpSx2* Dt* Evisco%specel(n)%visM* DumpSx3YZ                     
          else !zz

              ! MODIFIED
              LocStress(:,xx+ 3) =  Tdomain%specel(n)%DumpSx1* LocStress(:,xx+ 3) + &
                                    Tdomain%specel(n)%DumpSx2* sig2*dt

          endif

        endif
                  
        if  (Tdomain%rheononli)   then
          ! VISCO-ELASTOPLASTICITY
          !xz
          if (xx == 1)&
          Tdomain%specel(n)%gamma2(:,xx+ 3) = DumpSx3 

          !yz
          if (xx == 2)&
          Tdomain%specel(n)%gamma2(:,xx+ 3) = DumpSx3YZ 

          !zz
          if (xx == 3)&
          Tdomain%specel(n)%gamma2(:,xx+3) = DumpSx3ZZ
        endif
        deallocate(zipto)                     
        
      else if (.not. Tdomain%rheovisco .AND. Tdomain%rheononli)   then
        ! ELASTOPLASTICITY
        Tdomain%specel(n)%gamma2(:,xx+ 3)    = s2
      endif 
    endif 
    
    deallocate(s1,s2,hT,Vxloc,elVxloc,els1,els2)            
    deallocate(ni,lambda,mu)
    deallocate(sig1,sig2)

    ! Actual strain
    Tdomain%specel(n)%gamma2(:,xx+ 3) = Tdomain%specel(n)%gamma1(:,xx+ 3) + Dt* Tdomain%specel(n)%gamma2(:,xx+ 3) 
    !
  enddo

  ! Strain increment
  Tdomain%specel(n)%deps = Tdomain%specel(n)%gamma2- Tdomain%specel(n)%gamma1

enddo 

end subroutine sigeps



!=======================================================================
! ajout Celine

subroutine MAT_VISLA_strain2(zipt,dt,wp,ws,de1,de2,Mp,Ms,sig)

  double precision, intent(inout) :: zipt(8)
  double precision, intent(in) :: dt
  double precision, intent(in) :: Mp, Ms
  double precision, intent(in) :: wp(8), ws(8)
  double precision, intent(in) :: de1,de2
  double precision, intent(out) :: sig

  integer :: i
  double precision :: tau(8),zipto

  zipto = 0d0
  tau   = (/1.72333d-3,1.80701d-3,5.38887d-3,1.99322d-2,&
        8.49833d-2,4.09335d-1,2.05951d0,13.2629d0/)

  do i=1,8
    zipt(i)  = zipt(i)* exp(-dt/tau(i)) &
              + (wp(i) * Mp - 2. * ws(i) * Ms)* (1d0-exp(-dt/tau(i)))*(de1 + de2) &
              + 2. * ws(i)* (1d0-exp(-dt/tau(i)))*de1 * Ms
    zipto = zipto + zipt(i)
  enddo

  sig = ( Mp - 2. * Ms)* (de1 + de2)  + 2. * de1 * Ms - zipto

end subroutine MAT_VISLA_strain2

