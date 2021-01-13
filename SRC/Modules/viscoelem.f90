module viscoelem

implicit none

type :: viscoel
real, dimension(:),   pointer :: visM
real, dimension(:),   pointer :: DumpSx3
real, dimension(:),   pointer :: visw
real, dimension(:),   pointer :: viswp
real, dimension(:),   pointer :: visMp
real, dimension(:),   pointer :: DumpSx3YZ
real, dimension(:),   pointer :: DumpSx3ZZ
real, dimension(:),   pointer :: DumpSx3VOL

real, dimension(:,:), pointer :: visziptYZ
real, dimension(:,:), pointer :: visziptZZ
real, dimension(:,:), pointer :: visziptVOL
real, dimension(:,:), pointer :: viszipt

real  :: Q
real  :: fr
real  :: Qp 
end type

end module viscoelem