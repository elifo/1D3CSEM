module subdomains

type SubDomain
logical :: PML
logical :: Eeffective

integer	:: satsol
integer :: ngll,wpml
integer :: nonsol 	! Non-linearity

real :: dx, mu, dt , rho, Q, fr, Qp
real :: Emodulus, Gmodulus, Kmodulus, Vp, Vs, Ni
real :: lambda2mu
real :: midStress1, midStress2

real, dimension (:),allocatable :: GLLc,GLLpol,GLLw
real, dimension (:,:), allocatable :: hprime,hTprime


end type Subdomain

end module subdomains
