module nonelement
  
implicit none

type :: nonelem

   logical, dimension(:) ,pointer :: eliflogic

   integer, dimension(:), pointer :: aktifsur


   real, dimension(:) ,   pointer :: dTau
   real, dimension(:) ,   pointer :: qStress1
   real, dimension(:) ,   pointer :: qStress2
   real, dimension(:) ,   pointer :: p0Stress
   real, dimension(:) ,   pointer :: Gm0
   real, dimension(:) ,   pointer :: S0
   real, dimension(:) ,   pointer :: Sb 
   real, dimension(:) ,   pointer :: S
   real, dimension(:) ,   pointer :: S0old
   real, dimension(:) ,   pointer :: Sold
   real, dimension(:) ,   pointer :: Wn 
   real, dimension(:) ,   pointer :: r0                       
   real, dimension(:) ,   pointer :: r1                       
   real, dimension(:) ,   pointer :: Ws1
   real, dimension(:) ,   pointer :: Ws2
   real, dimension(:) ,   pointer :: Gmod1
   real, dimension(:) ,   pointer :: Gmod2
   real, dimension(:) ,   pointer :: SigmaMEff1
   real, dimension(:) ,   pointer :: SigmaMEff2 
   real, dimension(:) ,   pointer :: Emod1
   real, dimension(:) ,   pointer :: Emod2
   real, dimension(:) ,   pointer :: Kmod1
   real, dimension(:) ,   pointer :: Kmod2   
   real, dimension(:),    pointer :: Pore1
   real, dimension(:),    pointer :: Pore2
   real, dimension(:),    pointer :: gammaref 
   real, dimension(:),    pointer :: tauOCTO
   real, dimension(:),    pointer :: gamOCTO

   real, dimension(:,:),  pointer :: alf
   real, dimension(:,:),  pointer :: pbar
   real, dimension(:,:),  pointer :: hexaSigma
   real, dimension(:,:),  pointer :: hexaS2
   real, dimension(:,:),  pointer :: hexaS1
   real, dimension(:,:),  pointer :: dhexaS
   real, dimension(:,:),  pointer :: dhexaSigma
   real, dimension(:,:),  pointer :: plastF1
   real, dimension(:,:),  pointer :: plastF2
   real, dimension(:,:),  pointer :: dplastF1
   real, dimension(:,:),  pointer :: dplastF2                 
   real, dimension(:,:),  pointer :: SigmaEff1
   real, dimension(:,:),  pointer :: SigmaEff2
   real, dimension(:,:),  pointer :: SigmaP
   real, dimension(:,:),  pointer :: gammaSAT
   real, dimension(:,:),  pointer :: RSAT
   real, dimension(:,:),  pointer :: CNinvSAT

   real, dimension(:,:,:),pointer :: plastSa1
   real, dimension(:,:,:),pointer :: plastSa2
   real, dimension(:,:,:),pointer :: hexaE     

   real*8, dimension(:,:,:),pointer  :: hexaEd

end type

end module nonelement