subroutine define_mesh (Tdomain,meshfile,matfile)
  
use sdomain

implicit none

type (domain), intent(INOUT), target :: Tdomain
character (len=20), intent (IN)      :: matfile 
character (len=20), intent (IN)      :: meshfile


integer :: i
integer :: j
integer :: n
integer :: icount
integer :: nboundnode
integer :: nbmat
integer :: nummat
integer :: npml

integer, pointer  :: nel 
integer, pointer  :: ngllx
integer, pointer  :: ngllxb
integer, pointer  :: nnode
integer , pointer  :: ndomain        

character (len =1):: dummychar


! Reading meshfile

open (12,file=meshfile,status="old")
read (12,*)
read (12,*) 
read (12,*)
read (12,*)
read (12,*) ! Number of elements, Number of nodes, Number of domains

read (12,*)   Tdomain%nelem, Tdomain%nnode, Tdomain%nsubdomain
read (12,*) ! Nodes number and coordinates

nnode => Tdomain%nnode 
allocate (Tdomain%node(0:nnode-1))

Tdomain%borehole = .False.

nboundnode = 0
do i= 0,nnode-1
  read (12,*) j, Tdomain%node(j-1)%coord, Tdomain%node(j-1)%type 
  Tdomain%node(j-1)%num = j 
  if (Tdomain%node(j-1)%type == "D" .or. Tdomain%node(i)%type == "D") nboundnode = nboundnode+1

! BOREHOLE MODIFICATION
  if (Tdomain%node(j-1)%type == "B" .or. Tdomain%node(i)%type == "B") Tdomain%borehole = .True. 
enddo


if (Tdomain%borehole)  allocate(Tdomain%Vborehole(3))
! BOREHOLE MODIFICATION !



Tdomain%nbdirichlet = nboundnode
if (Tdomain%nbdirichlet > 0) then
  allocate(Tdomain%dirichlet(0:Tdomain%nbdirichlet-1))
  if (Tdomain%nbdirichlet == 1) then
    if (Tdomain%node(0)%type == "D") then
      Tdomain%dirichlet(0) = "L"
    else
      Tdomain%dirichlet(0) = "R"
    endif
  else
      Tdomain%dirichlet(0) = "L"
      Tdomain%dirichlet(1) = "R"
  endif
endif

read (12,*) ! Elements number, corresponding domain and nodes number
nel => Tdomain%nelem
ndomain => Tdomain%nsubdomain
allocate(Tdomain%specel(0:nel-1))

do i=0,nel-1
    allocate(Tdomain%specel(i)%Enode(0:1))
    read (12,*) j,Tdomain%specel(j-1)%Numdomain,&
                       Tdomain%specel(j-1)%Enode(0), Tdomain%specel(j-1)%Enode(1) 
    Tdomain%specel(j-1)%num = j
enddo

read (12,*) ! Domains number, corresponding  material number, dxmin and dxmax

allocate(Tdomain%dxmaxima(0:1,0:ndomain-1))
allocate(Tdomain%typemat (0:ndomain-1))
allocate(Tdomain%nbelem  (0:ndomain-1))

do i =0,ndomain-1
  read (12,*) j, Tdomain%typemat(j-1), Tdomain%nbelem(j-1), Tdomain%dxmaxima(0,j-1), Tdomain%dxmaxima(1,j-1)
enddo


! Boundary conditions

allocate(Tdomain%paramdiric(0:1,0:2))
if (nboundnode > 0) then
  read (12,*) ! Parameters for the Dirichlet conditions Left/Right (vmax,t0,tmax)
  do i=0,nboundnode-1
     read (12,*) dummychar,Tdomain%paramdiric(i,0),Tdomain%paramdiric(i,1),Tdomain%paramdiric(i,2)
  enddo
endif
close(12) 


! Reading material file

allocate(Tdomain%Sub_domain(0:Tdomain%nsubdomain-1))
open (13,file=matfile,status="unknown")
ndomain => Tdomain%nsubdomain
read (13,*)
read (13,*) ! Number of material types
read (13,*)   nbmat
read (13,*) ! Type Number, SVelocity, Density, PML, Delta T, Ngll, Q, Nonlinear Soil Type
npml = 0
do i=0,nbmat-1
  Tdomain%sub_domain(i)%wpml = -1
  read (13,*) j, Tdomain%Sub_domain(j-1)%Vp,Tdomain%Sub_domain(j-1)%Vs,Tdomain%Sub_domain(j-1)%rho, &
              Tdomain%Sub_domain(j-1)%PML, Tdomain%Sub_domain(j-1)%dt, Tdomain%Sub_domain(j-1)%ngll, &
              Tdomain%Sub_domain(j-1)%Qp, Tdomain%Sub_domain(j-1)%Q, &
              Tdomain%Sub_domain(j-1)%fr, Tdomain%Sub_domain(j-1)%nonsol
  
  ! Lame constants
  Tdomain%Sub_domain(j-1)%Mu        = Tdomain%Sub_domain(j-1)%rho* Tdomain%Sub_domain(j-1)%Vs**2 
  Tdomain%Sub_domain(j-1)%lambda2mu = Tdomain%Sub_domain(j-1)%rho* Tdomain%Sub_domain(j-1)%Vp**2

  if (Tdomain%Sub_domain(j-1)%PML) then
    Tdomain%sub_domain(j-1)%wpml = npml
    npml = npml + 1
  endif
enddo


Tdomain%n_pml = npml
if (npml > 0 ) then
  read(13,*); read(13,*)
  allocate(Tdomain%PMLCondition(0:npml-1))
  do i = 0,npml - 1
     read (13,*)  Tdomain%PMLCondition(i)%n_pow, Tdomain%PMLCondition(i)%A_Pow, Tdomain%PMLCondition(i)%Left
  
     Tdomain%PMLCondition(i)%Filtering = .False.
     Tdomain%PMLCondition(i)%freq      = 50
  enddo
else
  read(13,*); read(13,*); read(13,*)
endif

! Saturated/ Dry domains

if (Tdomain%Epressmod)  then
  read(13,*) ; read(13,*)
  j = 0
  do i = 0,nbmat-1
    Tdomain%Sub_domain(i)%Eeffective = .False.
    read(13,*)  j, Tdomain%Sub_domain(j- 1)%Eeffective, Tdomain%Sub_domain(j-1)%satsol  
  enddo
endif
close(13)

Tdomain%dtmaximum = 0
do i = 0,nbmat-1
    if (Tdomain%Sub_domain(i)%dt > Tdomain%dtmaximum ) & 
    Tdomain%dtmaximum = Tdomain%Sub_domain(i)%dt 
enddo


! Spread properties onto the elements
nel => Tdomain%nelem
do n = 0,nel-1

  nummat = Tdomain%specel(n)%Numdomain-1
  Tdomain%specel(n)%ngllx = Tdomain%Sub_domain(nummat)%ngll

  ngllx  => Tdomain%specel(n)%ngllx
  Tdomain%specel(n)%dt    = Tdomain%Sub_domain(nummat)%dt
  Tdomain%specel(n)%PML   = Tdomain%Sub_domain(nummat)%PML
  !
  Tdomain%specel(n)%Q     = Tdomain%Sub_domain(nummat)%Q
  Tdomain%specel(n)%Qp    = Tdomain%Sub_domain(nummat)%Qp
  Tdomain%specel(n)%fr    = Tdomain%Sub_domain(nummat)%fr



  if (Tdomain%specel(n)%PML) then
     npml = Tdomain%Sub_domain(nummat)%wpml
     Tdomain%specel(n)%FPML = .False.
     if (Tdomain%PMLCondition(npml)%Filtering) Tdomain%specel(n)%FPML = .True.
     Tdomain%specel(n)%Filter = exp (-Tdomain%PMLCondition(npml)%Freq * Tdomain%specel(n)%dt)
  endif
  allocate(Tdomain%specel(n)%Density  (0:ngllx-1))
  allocate(Tdomain%specel(n)%Mu       (0:ngllx-1))
  allocate(Tdomain%specel(n)%lambda2mu(0:ngllx-1))

  nummat = Tdomain%specel(n)%Numdomain-1
  do i = 0, ngllx-1
      Tdomain%specel(n)%Density(i)   = Tdomain%Sub_domain(nummat)%rho
      Tdomain%specel(n)%Mu(i)        = Tdomain%Sub_domain(nummat)%Mu
      Tdomain%specel(n)%lambda2mu(i) = Tdomain%Sub_domain(nummat)%lambda2mu
  enddo
enddo


nel => Tdomain%nelem
icount = 0

do n=0,nel-1
! Global numerotation for GLL points
  ngllx  => Tdomain%specel(n)%ngllx
  allocate(Tdomain%specel(n)%Iglobnum(0:ngllx-1)) 

  do i = 0,ngllx-1
    Tdomain%specel(n)%Iglobnum(i) = icount 
    icount = icount+1

    if (i == 0 .and. n>0) then
        ngllxb => Tdomain%specel(n-1)%ngllx
        Tdomain%specel(n)%Iglobnum(i) = Tdomain%specel(n-1)%Iglobnum(ngllxb-1)
        icount = icount- 1
    endif
  enddo
enddo
Tdomain%npoints = icount


! Define elements connected at each node

do i = 0,nnode-1
  allocate(Tdomain%node(i)%elemconnect(0:1))
  Tdomain%node(i)%elemconnect = -1
  do n = 0,nel-1
    if (Tdomain%specel(n)%Enode(0) == Tdomain%node(i)%num) Tdomain%node(i)%elemconnect(1) = n 
    if (Tdomain%specel(n)%Enode(1) == Tdomain%node(i)%num) Tdomain%node(i)%elemconnect(0) = n 
  enddo
enddo

return
end subroutine define_mesh
