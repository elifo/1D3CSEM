module stimeparam 

! This Module is a block with time parameters

type :: time
integer :: nt,skips,Newmark,itermax,output
real :: beta, gamm, dt,rtime                          
real :: alpha, duration
end type

end module stimeparam
