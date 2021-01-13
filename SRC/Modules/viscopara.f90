module viscopara

use viscoelem

implicit none

type :: visco
real, dimension(8) :: visal
real, dimension(8) :: visbe 
real, dimension(8) :: vistau

type(viscoel), dimension(:), pointer :: specel
end type

contains

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine viscoe(Evisco,Tdomain)

use sdomain

implicit none

type(visco),  intent(INOUT), target :: Evisco
type(domain), intent(IN),    target :: Tdomain


complex     :: A
complex     :: Ap 

integer     :: i 
integer     :: n
integer     :: nummat

real        :: chi
real        :: chip

integer, pointer      :: nel 
integer(kind=8), pointer      :: ngllx
integer, pointer      :: ndom

complex, dimension(8) :: A1
complex, dimension(8) :: A2 
complex, dimension(8) :: A1p
complex, dimension(8) :: A2p

!
! [Liu & Archuleta, 2006 - Eqns 1-8 & Table 1]


Evisco%visal  = 0.0
Evisco%visbe  = 0.0
Evisco%vistau = 0.0 

Evisco%visal  = (/1.66958e-2,3.81644e-2,9.84666e-3,-1.36803e-2,-2.85125e-2,-5.37309e-2,-6.65035e-2,-1.33696e-1/)
Evisco%visbe  = (/8.98758e-2,6.84635e-2,9.67052e-2,1.20172e-1,1.30728e-1,1.38746e-1,1.40705e-1,2.14647e-1/)
Evisco%vistau = (/1.72333e-3,1.80701e-3,5.38887e-3,1.99322e-2,8.49833e-2,4.09335e-1,2.05951,13.2629/)


nel   => Tdomain%nelem
allocate(Evisco%specel(0:nel-1)) 

! Memory variables calculation for shear and axial components

do n = 0,nel-1
  allocate(Evisco%specel(n)%visw (8))
  allocate(Evisco%specel(n)%viswp(8))

  Evisco%specel(n)%Qp  = Tdomain%specel(n)%Qp 
  Evisco%specel(n)%Q   = Tdomain%specel(n)%Q 
  Evisco%specel(n)%fr  = Tdomain%specel(n)%fr 

  chi  = (3.071+ 1.433*(Evisco%specel(n)%Q** (-1.158))*log(real(Evisco%specel(n)%Q/ 5.0)))/(1.0+ 0.415*Evisco%specel(n)%Q)
  chip = (3.071+ 1.433*(Evisco%specel(n)%Qp**(-1.158))*log(real(Evisco%specel(n)%Qp/5.0)))/(1.0+ 0.415*Evisco%specel(n)%Qp)

  do i = 1,8
    Evisco%specel(n)%visw(i)  = chi*  (chi*  Evisco%visal(i)+ Evisco%visbe(i))
    Evisco%specel(n)%viswp(i) = chip* (chip* Evisco%visal(i)+ Evisco%visbe(i))
  end do

end do


print*,"Viscoelastic Modulus Calculation"
! Viscoelastic unrelaxed modulus Mu: shear components

do n = 0,nel-1

  do i = 1,8
    A1(i) = 1.0+ cmplx(0.0,1.0)*2.0*(4.0*atan(1.0))*Evisco%vistau(i)* Evisco%specel(n)%fr  
    A2(i) = Evisco%specel(n)%visw(i)/A1(i)
  end do

  A = A2(1)
  do i = 2,8
    A = A+ A2(i)
  end do

  ngllx => Tdomain%specel(n)%ngllx
  allocate(Evisco%specel(n)%visM(0:ngllx-1))
  do i = 0,ngllx-1
    Evisco%specel(n)%visM(i) = Tdomain%specel(n)%Mu(i)/cabs(1.0-A)
  end do

end do


! Viscoelastic unrelaxed modulus Mu: axial component

do n = 0,nel-1
  do i = 1,8
    A1p(i) = 1.0+ cmplx(0.0,1.0)*2.0*(4.0*atan(1.0))*Evisco%vistau(i)* Evisco%specel(n)%fr  
    A2p(i) = Evisco%specel(n)%viswp(i)/A1p(i)
  end do

  Ap = A2p(1)
  do i = 2,8
    Ap = Ap+ A2p(i)
  end do

  ngllx => Tdomain%specel(n)%ngllx
  allocate(Evisco%specel(n)%visMp(0:ngllx-1))
  do i = 0,ngllx-1
    Evisco%specel(n)%visMp(i) = Tdomain%specel(n)%lambda2mu(i)/cabs(1.0-Ap)
  end do
end do

end subroutine viscoe
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module viscopara
