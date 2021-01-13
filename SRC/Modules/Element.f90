module selement
! This module contains the derived type: element, which contains the 
! coordinates of the Gauss-Lobatto-Legendre, the Displacement, Velocity and
! Axeleration in each GLL node, the elastic properties, the external forces,
! the mass matrix, the internal forces, the Jacobian, the inverse Gradient,
! the internal coefficients and the coordinates of the control points.
! 
implicit none


type :: element
integer :: num,ngllx,Numdomain
integer, dimension(:), pointer :: Enode
integer, dimension (:), pointer :: Iglobnum
real :: dt
!real, dimension(:), pointer ::      
!real, dimension(:), pointer ::     
real, dimension (:), pointer :: Density, Mu, lambda2mu
real, dimension(:), pointer :: MassMat, Jacob
real, dimension(:), pointer :: InvGrad
real, dimension(:), pointer :: ACoeff
real :: Q, fr, Qp
real, dimension(:,:), pointer :: Stress, gamma1, gamma2, Veloc, V0,Calc0,Accel,Forces,MidDispl
real, dimension(:,:), pointer :: SigmaIni


! 1D3C
real :: Emodulus, Gmodulus, Kmodulus, Vp, Vs, Ni
real, dimension(:), pointer :: epsm, sigm, depsm,dsigm
real, dimension(:,:), pointer :: deps, dsigma, de, dplastic


!!!Incident Wave
real, dimension(:), pointer :: elStress, elMidDispl


! Fault requiring
!integer :: nfault

! PML allocation 
logical :: PML,FPML
real :: Filter
real, dimension (:), pointer :: DumpMass
real, dimension (:), pointer :: DumpSx1,DumpSx2,DumpVx1,DumpVx2
real, dimension (:), pointer :: Is, Iv, Istress, Iveloc  
real, dimension (:), pointer :: DumpMassE
real, dimension (:), pointer :: DumpSx1E,DumpSx2E,DumpVx1E,DumpVx2E   
real, dimension (:), pointer :: IsE, IvE

end type

end module selement
