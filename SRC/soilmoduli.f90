subroutine soilmoduli(Tdomain)

use sdomain

implicit none

type (domain), intent(INOUT), target :: Tdomain

integer           ::  i
integer           ::  nummat
integer           ::  n

integer(kind=8), pointer  ::  numdom 

real, pointer     ::  Gmodulus
real, pointer     ::  Ni
real, pointer     ::  Emodulus
real, pointer     ::  Kmodulus
real, pointer     ::  Vp
real, pointer     ::  Vs
real, pointer     ::  Mu


numdom    =>    Tdomain%nsubdomain    

do i = 0,numdom- 1

  Gmodulus      =>    Tdomain%Sub_domain(i)%Gmodulus 
  Ni            =>    Tdomain%Sub_domain(i)%Ni
  Emodulus      =>    Tdomain%Sub_domain(i)%Emodulus
  Kmodulus      =>    Tdomain%Sub_domain(i)%Kmodulus 
  Vp            =>    Tdomain%Sub_domain(i)%Vp
  Vs            =>    Tdomain%Sub_domain(i)%Vs
  Mu            =>    Tdomain%Sub_domain(i)%Mu 

  if ((Vp/Vs) == 1.0)   &
  print*,   "POISSON RATIO = Infinity ! Check P and S Velocities" 

  Gmodulus      =   Mu 
  Ni            =   (Vp/ Vs)**2.0
  Ni            =   (Ni- 2.0)/ (2.0*(Ni- 1.0))

  Emodulus      =   2.0* Gmodulus* (1.0+ Ni)
  Kmodulus      =   Emodulus/ (3.0*(1.0- 2.0* Ni)) 
enddo 


! Spread properties onto the elements

do n = 0,Tdomain%nelem-1
  nummat        =   Tdomain%specel(n)%Numdomain- 1

  Tdomain%specel(n)%Gmodulus  = Tdomain%Sub_domain(nummat)%Gmodulus
  Tdomain%specel(n)%Emodulus  = Tdomain%Sub_domain(nummat)%Emodulus
  Tdomain%specel(n)%Kmodulus  = Tdomain%Sub_domain(nummat)%Kmodulus
  Tdomain%specel(n)%Ni        = Tdomain%Sub_domain(nummat)%Ni            
enddo

end subroutine soilmoduli