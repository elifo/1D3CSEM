subroutine compute_GLL(Tdomain)

use sdomain

implicit none
type (Domain), intent (INOUT) :: Tdomain

! Local declarations 
integer :: i,ndomains
integer, pointer ::  ngll

ndomains = Tdomain%nsubdomain

do i = 0, ndomains-1
  ngll => Tdomain%Sub_Domain(i)%ngll

  allocate (Tdomain%Sub_Domain(i)%GLLc (0:ngll-1))
  allocate (Tdomain%Sub_Domain(i)%GLLpol (0:ngll-1))
  allocate (Tdomain%Sub_Domain(i)%GLLw (0:ngll-1))
  allocate (Tdomain%Sub_Domain(i)%hprime (0:ngll-1,0:ngll-1))
  allocate (Tdomain%Sub_Domain(i)%hTprime (0:ngll-1,0:ngll-1))


 ! USING FUNARO SUBROUTINES
 ! ZELEGL computes the coordinates of GLL points
 ! WELEGL computes the respective weights
 ! DMLEGL compute the matrix of the first derivatives in GLL points
  call zelegl (ngll-1,Tdomain%Sub_Domain(i)%GLLc,Tdomain%Sub_Domain(i)%GLLpol)
  call welegl (ngll-1, Tdomain%Sub_Domain(i)%GLLc, Tdomain%Sub_Domain(i)%GLLpol, Tdomain%Sub_Domain(i)%GLLw)


  call dmlegl (ngll-1, ngll-1, Tdomain%Sub_Domain(i)%GLLc, Tdomain%Sub_Domain(i)%GLLpol, Tdomain%Sub_Domain(i)%hTprime)
  Tdomain%Sub_Domain(i)%hprime =  TRANSPOSE ( Tdomain%Sub_Domain(i)%hTprime )

  deallocate(Tdomain%Sub_Domain(i)%GLLpol)
enddo
return
end subroutine compute_GLL
