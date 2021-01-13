subroutine Define_Arrays(Tdomain)

use sdomain
use logical_input

implicit none

type (domain),intent (INOUT), target :: Tdomain


integer  :: i
integer  :: n
integer  :: n1
integer  :: n2
integer  :: i_domain
integer  :: npow
integer  :: n_pml
integer  :: elif=1

real     :: c0
real     :: vs
real     :: ri
real     :: dx 
real     :: Apow
real     :: vp
real, external  :: pow

integer, pointer     :: nel
integer, pointer     :: ngllx

real, dimension (:), pointer        :: Jac
real, dimension (:), pointer        :: Rmu
real, dimension (:), pointer        :: Rlam

real, dimension (:), allocatable    :: Whei,Id  
real, dimension (:), allocatable    :: wx
real, dimension (:), allocatable    :: wx2




nel => Tdomain%nelem

do n = 0,nel-1

  ngllx => Tdomain%specel(n)%ngllx
 
  Jac  => Tdomain%specel(n)%Jacob
  Rmu  => Tdomain%specel(n)%Mu
  Rlam => Tdomain%specel(n)%lambda2mu
  

!
  allocate(Whei (0:ngllx-1))
  allocate(wx   (0:ngllx-1))
  allocate(Id   (0:ngllx-1))
  allocate(wx2  (0:ngllx-1))
!
  Id = 1.0
  i_domain = Tdomain%specel(n)%NumDomain-1

  Whei     = Tdomain%Sub_Domain(i_domain)%GLLw

  Tdomain%specel(n)%MassMat = Whei*Tdomain%specel(n)%Density* Jac
  Tdomain%specel(n)%Acoeff  = -Tdomain%specel(n)%InvGrad* Whei * Jac
 !
  ! Calculation of wx
  if (.not. Tdomain%specel(n)%PML) then
    wx  = 0.0
    wx2 = 0.0
  else
    dx    = Tdomain%node(Tdomain%specel(n)%Enode(1)-1)%coord - &
            Tdomain%node(Tdomain%specel(n)%Enode(0)-1)%coord
    n_pml = Tdomain%Sub_domain(i_domain)%wpml
    Apow  = Tdomain%PMLCondition(n_pml)%A_pow
    npow  = Tdomain%PMLCondition(n_pml)%n_pow

    if (Tdomain%specel(n)%FPML) c0 = Tdomain%specel(n)%Filter
    if (.not. Tdomain%PMLCondition(n_pml)%Left) then
      do i = 0,ngllx-1
        ri = 0.5*(1+Tdomain%Sub_Domain(i_domain)%GLLc(i))*float (ngllx-1)
        ! Shear components
        vs    = Rmu(i)/Tdomain%specel(n)%Density(i)
        vs    = sqrt(vs)
        wx(i) = pow(ri,vs,ngllx-1,dx,npow,Apow)
       
        ! axial component
        vp    = Rlam(i)/ Tdomain%specel(n)%Density(i)
        vp    = sqrt(vp)            
        wx2(i)= pow(ri,vp,ngllx-1,dx,npow,Apow)
      enddo
    else
      do i = 0,ngllx-1
        ri = 0.5*(1+Tdomain%Sub_Domain(i_domain)%GLLc(ngllx-1-i))* float (ngllx-1)
        ! Shear components
        vs =  Rmu(i)/Tdomain%specel(n)%Density(i)
        vs = sqrt(vs)
        wx(i) = pow (ri,vs,ngllx-1,dx,npow,Apow)
       
        ! axial component
        vp =  Rlam(i)/Tdomain%specel(n)%Density(i)
        vp = sqrt(vp)            
        wx2(i) = pow (ri,vp,ngllx-1,dx,npow,Apow)
      enddo
    endif    
  end if

! Strong formulation for stress
! Shear components
  if (Tdomain%specel(n)%FPML) then
    Tdomain%specel(n)%DumpSx2  = Id + 0.5 * Tdomain%specel(n)%dt * wx * c0
    Tdomain%specel(n)%DumpSx1  = Id - 0.5 * Tdomain%specel(n)%dt * wx *  c0
    Tdomain%specel(n)%DumpSx2  = 1.0/ Tdomain%specel(n)%DumpSx2
    Tdomain%specel(n)%DumpSx1  = Tdomain%specel(n)%DumpSx1 * Tdomain%specel(n)%DumpSx2
    Tdomain%specel(n)%Is       =  Tdomain%specel(n)%dt * wx * c0
    Tdomain%specel(n)%DumpVx2  = 0.5 *Tdomain%specel(n)%dt * Tdomain%specel(n)%Density * Whei * wx * Jac * c0
    Tdomain%specel(n)%DumpVx1  = 0.5 *Tdomain%specel(n)%dt * Tdomain%specel(n)%Density * Whei * wx * Jac * c0
    Tdomain%specel(n)%IV       = Tdomain%specel(n)%Density * Whei * Jac * Tdomain%specel(n)%dt * wx * c0
    Tdomain%specel(n)%DumpMass = Tdomain%specel(n)%Density * Whei * wx * Jac * c0
  else 
    Tdomain%specel(n)%DumpSx2  = Id + 0.5 * Tdomain%specel(n)%dt * wx 
    Tdomain%specel(n)%DumpSx2  = 1.0/ Tdomain%specel(n)%DumpSx2
    Tdomain%specel(n)%DumpSx1  = (Id - 0.5 * Tdomain%specel(n)%dt * wx) * Tdomain%specel(n)%DumpSx2
    Tdomain%specel(n)%DumpMass = Tdomain%specel(n)%Density * Whei * wx * Jac
  endif

! Axial component
  if (Tdomain%specel(n)%FPML) then
    Tdomain%specel(n)%DumpSx2E  = Id + 0.5 * Tdomain%specel(n)%dt * wx2 * c0
    Tdomain%specel(n)%DumpSx1E  = Id - 0.5 * Tdomain%specel(n)%dt * wx2 * c0
    Tdomain%specel(n)%DumpSx2E  = 1.0/ Tdomain%specel(n)%DumpSx2E
    Tdomain%specel(n)%DumpSx1E  = Tdomain%specel(n)%DumpSx1E * Tdomain%specel(n)%DumpSx2E
    Tdomain%specel(n)%IsE       = Tdomain%specel(n)%dt * wx2 * c0
    Tdomain%specel(n)%DumpVx2E  = 0.5 *Tdomain%specel(n)%dt * Tdomain%specel(n)%Density * Whei * wx2 * Jac * c0
    Tdomain%specel(n)%DumpVx1E  = 0.5 *Tdomain%specel(n)%dt * Tdomain%specel(n)%Density * Whei * wx2 * Jac * c0
    Tdomain%specel(n)%IVE       = Tdomain%specel(n)%Density * Whei * Jac * Tdomain%specel(n)%dt * wx2 * c0
    Tdomain%specel(n)%DumpMassE = Tdomain%specel(n)%Density * Whei * wx2 * Jac * c0
  else 
    Tdomain%specel(n)%DumpSx2E  = Id + 0.5 * Tdomain%specel(n)%dt * wx2 
    Tdomain%specel(n)%DumpSx2E  = 1.0/ Tdomain%specel(n)%DumpSx2E
    Tdomain%specel(n)%DumpSx1E  = (Id - 0.5 * Tdomain%specel(n)%dt * wx2) * Tdomain%specel(n)%DumpSx2E
    Tdomain%specel(n)%DumpMassE = Tdomain%specel(n)%Density * Whei * wx2 * Jac
  endif

  deallocate(Whei,wx,Id) 
  deallocate(wx2)
enddo

! Communicating mass matrix between the coefficients

do i = 0,Tdomain%nnode-1
  if (Tdomain%node(i)%type == "I") then 
    n1 = Tdomain%node(i)%elemconnect(0)
    n2 = Tdomain%node(i)%elemconnect(1)
    
    ngllx => Tdomain%specel(n1)%ngllx
    
    Tdomain%specel(n2)%MassMat(0)       = Tdomain%specel(n2)%MassMat(0)+Tdomain%specel(n1)%MassMat(ngllx-1)
    Tdomain%specel(n1)%MassMat(ngllx-1) = Tdomain%specel(n2)%MassMat(0)
    if (Tdomain%specel(n1)%FPML .and. Tdomain%specel(n2)%FPML) then
      Tdomain%specel(n2)%DumpMass(0)       = Tdomain%specel(n2)%DumpMass(0)+Tdomain%specel(n1)%DumpMass(ngllx-1)
      Tdomain%specel(n1)%DumpMass(ngllx-1) = Tdomain%specel(n2)%DumpMass(0)
      Tdomain%specel(n2)%DumpVx1(0)        = Tdomain%specel(n2)%DumpVx1(0)+Tdomain%specel(n1)%DumpVx1(ngllx-1)
      Tdomain%specel(n1)%DumpVx1(ngllx-1)  = Tdomain%specel(n2)%DumpVx1(0)
      Tdomain%specel(n2)%DumpVx2(0)        = Tdomain%specel(n2)%DumpVx2(0)+Tdomain%specel(n1)%DumpVx2(ngllx-1)
      Tdomain%specel(n1)%DumpVx2(ngllx-1)  = Tdomain%specel(n2)%DumpVx2(0)
      Tdomain%specel(n2)%IV(0)             = Tdomain%specel(n2)%IV(0)+Tdomain%specel(n1)%IV(ngllx-1)
      Tdomain%specel(n1)%IV(ngllx-1)       = Tdomain%specel(n2)%IV(0)

      Tdomain%specel(n2)%DumpMassE(0)      = Tdomain%specel(n2)%DumpMassE(0)+Tdomain%specel(n1)%DumpMassE(ngllx-1)
      Tdomain%specel(n1)%DumpMassE(ngllx-1)= Tdomain%specel(n2)%DumpMassE(0)
      Tdomain%specel(n2)%DumpVx1E(0)       = Tdomain%specel(n2)%DumpVx1E(0)+Tdomain%specel(n1)%DumpVx1E(ngllx-1)
      Tdomain%specel(n1)%DumpVx1E(ngllx-1) = Tdomain%specel(n2)%DumpVx1E(0)
      Tdomain%specel(n2)%DumpVx2E(0)       = Tdomain%specel(n2)%DumpVx2E(0)+Tdomain%specel(n1)%DumpVx2E(ngllx-1)
      Tdomain%specel(n1)%DumpVx2E(ngllx-1) = Tdomain%specel(n2)%DumpVx2E(0)
      Tdomain%specel(n2)%IVE(0)            = Tdomain%specel(n2)%IVE(0)+Tdomain%specel(n1)%IVE(ngllx-1)
      Tdomain%specel(n1)%IVE(ngllx-1)      = Tdomain%specel(n2)%IVE(0)    
    endif
  endif
enddo


! Mat PML

do n = 0, nel - 1  
  if (Tdomain%specel(n)%FPML) then
    Tdomain%specel(n)%DumpVx2  = Tdomain%specel(n)%MassMat + Tdomain%specel(n)%DumpVx2
    Tdomain%specel(n)%DumpVx2  = 1.0/Tdomain%specel(n)%DumpVx2 
    Tdomain%specel(n)%DumpVx1  = (Tdomain%specel(n)%MassMat - Tdomain%specel(n)%DumpVx1) * Tdomain%specel(n)%DumpVx2

    Tdomain%specel(n)%DumpVx2E = Tdomain%specel(n)%MassMat + Tdomain%specel(n)%DumpVx2E
    Tdomain%specel(n)%DumpVx2E = 1.0/Tdomain%specel(n)%DumpVx2E
    Tdomain%specel(n)%DumpVx1E = (Tdomain%specel(n)%MassMat - Tdomain%specel(n)%DumpVx1E) * Tdomain%specel(n)%DumpVx2E
  else
    Tdomain%specel(n)%DumpVx2 = Tdomain%specel(n)%MassMat + 0.5 *Tdomain%specel(n)%dt * Tdomain%specel(n)%DumpMass 
    Tdomain%specel(n)%DumpVx2 = 1.0/Tdomain%specel(n)%DumpVx2 
    Tdomain%specel(n)%DumpVx1 = (Tdomain%specel(n)%MassMat - 0.5 * &
        Tdomain%specel(n)%dt * Tdomain%specel(n)%DumpMass) * Tdomain%specel(n)%DumpVx2

!!!
    Tdomain%specel(n)%DumpVx2E = Tdomain%specel(n)%MassMat + 0.5 *Tdomain%specel(n)%dt * Tdomain%specel(n)%DumpMassE 
    Tdomain%specel(n)%DumpVx2E = 1.0/Tdomain%specel(n)%DumpVx2E
    Tdomain%specel(n)%DumpVx1E = (Tdomain%specel(n)%MassMat - 0.5 * &
        Tdomain%specel(n)%dt * Tdomain%specel(n)%DumpMassE) * Tdomain%specel(n)%DumpVx2E         
  endif     
enddo

return

end subroutine
