module sreceiver_properties

type :: receiver

integer 	:: SSnr
integer 	:: Nr 

real        :: xRec
real        :: xi
real        :: SSxRec
real        :: SSxi

real, dimension(:), pointer :: Interp_coeff
real, dimension(:), pointer :: SSInterp_coeff

end type 

end module sreceiver_properties
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       

module sreceivers

use sreceiver_properties

type :: receiver_stream

integer 	:: nsta
integer 	:: nsta1
integer 	:: SSnsta

real, dimension(:,:), pointer :: StoreTraceAzz
real, dimension(:,:), pointer :: StoreTraceVzz
real, dimension(:,:), pointer :: SSStoreTrace1zz
real, dimension(:,:), pointer :: SSStoreTrace2zz

real, dimension(:,:), pointer :: StoreTraceAxz
real, dimension(:,:), pointer :: StoreTraceVxz
real, dimension(:,:), pointer :: SSStoreTrace1xz
real, dimension(:,:), pointer :: SSStoreTrace2xz

real, dimension(:,:), pointer :: StoreTraceAyz
real, dimension(:,:), pointer :: StoreTraceVyz
real, dimension(:,:), pointer :: SSStoreTrace1yz
real, dimension(:,:), pointer :: SSStoreTrace2yz

real, dimension(:,:), pointer :: StoreS
real, dimension(:,:), pointer :: StoreS0
real, dimension(:,:), pointer :: StoreT
real, dimension(:,:), pointer :: Storer
real, dimension(:,:), pointer :: StoreWs
real, dimension(:,:), pointer :: StoreWn
real, dimension(:,:), pointer :: StorePore
real, dimension(:,:), pointer :: StoreGref
real, dimension(:,:), pointer :: StoreSOcto
real, dimension(:,:), pointer :: StoreGOcto
real, dimension(:,:), pointer :: StoreGmod
real, dimension(:,:), pointer :: StoreEmod
real, dimension(:,:), pointer :: StoreKmod
real, dimension(:,:), pointer :: StoreGm0
real, dimension(:,:), pointer :: StoreNi
real, dimension(:,:), pointer :: StoreGact


type (receiver) , dimension(:), pointer :: Rec
type (receiver) , dimension(:), pointer :: SSRec

end type 

end module sreceivers