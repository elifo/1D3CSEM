subroutine pol_lagrange (n,GLLc,k,x,y)
implicit none
integer, intent (IN) :: n,k
real, dimension  (0:n-1), intent (IN) :: GLLc
real,intent (IN) :: x
real, intent (OUT) :: y

integer :: i

y = 1
if (n ==0 ) then
   write(*,"(A)", ADVANCE = "NO") "Bad n number - Press <Return> to continue."
   read(*,"(A)") 
endif
if (n ==1 ) return
   do i =0,n-1
      if (i /= k)  y = y * (x-GLLc(i))/(GLLc(k)-GLLc(i))
   enddo
return
end subroutine pol_lagrange
  
