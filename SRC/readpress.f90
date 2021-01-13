subroutine readpress(Enonli)
   
use nonlinear

implicit none

type(nonli), intent(INOUT),target :: Enonli

integer :: numdom
integer :: i 

double precision :: phi, phi_p

! Reading saturation parameters

open(61, file = "pressparam.dat", status = "old")
read(61,*) 
read(61,*) numdom
    
allocate(Enonli%SATm1(numdom))
allocate(Enonli%SATm2(numdom))
allocate(Enonli%SATp1(numdom))
allocate(Enonli%SATp2(numdom))
allocate(Enonli%SATS1(numdom))
allocate(Enonli%SATw1(numdom))
allocate(Enonli%SATc1(numdom)) 
allocate(Enonli%Cohesion(numdom))
allocate(Enonli%Isoconso(numdom)) 
allocate(Enonli%cos_phif(numdom)) !!!


Enonli%SATm1 = 0.0
Enonli%SATm2 = 0.0
Enonli%SATp1 = 0.0
Enonli%SATp2 = 0.0
Enonli%SATS1 = 0.0
Enonli%SATw1 = 0.0
Enonli%SATc1 = 0.0 
Enonli%Cohesion = 0.0
Enonli%Isoconso = 0.0 
Enonli%cos_phif = 0.0


do i = 1,numdom
  read(61,*)
  read(61,*) phi
  phi = phi* 4.0*atan(1.0)/ 180.0
  Enonli%SATm1(i) = dsin(phi)
 
  !!!
  Enonli%cos_phif(i) = dcos(phi)

  read(61,*) phi_p
  phi_p = phi_p* 4.0*atan(1.0)/ 180.0
  Enonli%SATm2(i) = dsin(phi_p)


  read(61,*) Enonli%SATp1(i)
  read(61,*) Enonli%SATp2(i)
  read(61,*) Enonli%SATs1(i)
  read(61,*) Enonli%SATw1(i)
  read(61,*) Enonli%Cohesion(i)
  read(61,*) Enonli%Isoconso(i)
enddo  


end subroutine readpress