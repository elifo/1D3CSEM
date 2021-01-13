module nonlinear

use nonelement

type :: nonli
  logical :: modeltype

  integer :: Nsoil
  integer :: Nspr

  integer, dimension(:), pointer :: Npg
  integer, dimension(:), pointer :: soil_type

  real, dimension(:),   pointer  :: gammaref
  real, dimension(:),   pointer  :: SATm1
  real, dimension(:),   pointer  :: cos_phif  !!!

  real, dimension(:),   pointer  :: SATm2
  real, dimension(:),   pointer  :: SATp1
  real, dimension(:),   pointer  :: SATp2
  real, dimension(:),   pointer  :: SATS1
  real, dimension(:),   pointer  :: SATw1
  real, dimension(:),   pointer  :: SATc1
  real, dimension(:),   pointer  :: Cohesion
  real, dimension(:),   pointer  :: Isoconso


  real, dimension(:,:), pointer  :: straing
  real, dimension(:,:), pointer  :: GoverG
  real, dimension(:,:), pointer  :: H
  real, dimension(:,:), pointer  :: R
  real, dimension(:,:), pointer  :: CNinv

  type(nonelem), dimension(:), pointer :: specel
end type





contains
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine get_G_info (Enonli)

  implicit none

  type(nonli), intent(INOUT), target     ::  Enonli

  integer        ::    k 
  integer        ::    kg
  character*20   ::    fnamef

  ! Reading the soil G/Gmax data
  open(61, file="GoverGmax.dat", status ="old")
  read(61,*) 
  read(61,*)  Enonli%modeltype
  read(61,*)  Enonli%Nspr
  read(61,*)
  if ( Enonli%modeltype ) then		                         	! If hyperbolic   NL model :
    read(61,*)  Enonli%Nsoil
    allocate(Enonli%gammaref(Enonli%Nsoil))	
    do k=1,Enonli%Nsoil
	   	read(61,*)  Enonli%gammaref(k)
    enddo
  else                                            					! If experimental NL model :
!
    read(61,*)  Enonli%Nsoil
    do k=1,Enonli%Nsoil
      read(61,*)  
    enddo
    read(61,*) ; read(61,*)
!
    read(61,*)  Enonli%Nsoil
   	allocate(Enonli%straing(Enonli%Nspr,100))			
	  allocate(Enonli%GoverG(Enonli%Nspr,100))			
	  allocate(Enonli%Npg(Enonli%Nsoil))				
	  allocate(Enonli%soil_type(Enonli%Nsoil)) 			

  	do k = 1,Enonli%Nsoil
	    read(61,*)  Enonli%soil_type(k)
	    read(61,*)  Enonli%Npg(k)
	    do kg = 1,Enonli%Npg(k)
        read(61,*)  Enonli%straing(kg,k), Enonli%GoverG(kg,k)
        Enonli%straing(kg,k) = Enonli%straing(kg,k)/100.0	           
      enddo
	  enddo
  endif
  close(61) 
  
  !write(77,*) Enonli%modeltype,Enonli%Nspr,Enonli%Nsoil,Enonli%gammaref
  print*, "SPRING NUMBER   ",    Enonli%Nspr
  print*, "SOIL TYPE NUMBER",    Enonli%Nsoil

end subroutine get_G_info
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
subroutine nonline_curve(Tdomain,Enonli,Gmax,isoil,R,H,CNinv)
   
  use sdomain
   
  implicit none
  
  type(nonli),    intent(INOUT),  target   :: Enonli
  type(domain),   intent(IN),     target   :: Tdomain
  integer,        intent(IN)               :: isoil   
  real*8,           intent(IN)               :: Gmax
   
  real      :: x     (Enonli%Npg(isoil))		
  real      :: GoverG(Enonli%Npg(isoil)) 
  real      :: y2    (Enonli%Npg(isoil))
  real      :: strain(Enonli%Npg(isoil))
  integer   :: i,j,n, nspr			
  real      :: x0,xu,dx
  real      :: gamma (Enonli%Nspr)
  real      :: Gout,elo
  real      :: G(Enonli%Nspr) 
  real      :: R(Enonli%Nspr),H(Enonli%Nspr), CNinv(Enonli%Nspr- 1)

  nspr = Enonli%Nspr

  i = isoil
  do j = 1,Enonli%Npg(i)
    x(j)      = alog10(Enonli%straing(j,i))
    GoverG(j) = Enonli%GoverG(j,i)
    strain(j) = Enonli%straing(j,i)
  enddo
 
  n = Enonli%Npg(isoil)
  !call spline(x,GoverG,n,0.0,0.0,y2)

  x0 = log10(strain(1))
  xu = log10(strain(n))
  dx = (xu-x0)/(Enonli%Nspr- 1) 

  do i = 1,nspr
    gamma(i) = 10.0**(x0 + dx*(i- 1))
    !call splint(x,GoverG,y2,n,alog10(sngl(gamma(i))),Gout)
    call hyper_interp(strain,GoverG,n,sngl(gamma(i)),Gout)
    G(i) = Gout
    R(i) = Gmax* G(i)* gamma(i)
    write(33,*)   gamma(i), R(i)	                
  enddo

  call CNinverse(Gmax,nspr,gamma,R,CNinv)

  return
end subroutine nonline_curve

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine hyper_interp(x,y,n,xi,yi)
!	A PROGRAM FOR INTERPOLATION BY ARC HYPERBOLIC
!	SINE METHOD BETWEEN TWO POINTS
!	x--ABSCISSA VALUE OF INPUT POINTS;  y--ORDINATE
! 	VALUE OF INPUT POINTS;  xi--INPUT ABSCISSA VALUE;
!	yi--OUTPUT ORDINATE VALUE

   dimension x(n), y(n)

	 if (x(1)-1.e-10 < 1.e-10) x(1)=1.e-10
   if (xi <= x(1)) then
      yi=y(1)
      return
   else
      j=2
   end if
30 if (xi < x(j)) then
	    ars  = alog(xi*1.e+4+((xi*1.e+4)**2+1.)**0.5)
	    ars1 = alog(x(j-1)*1.e+4+((x(j-1)*1.e+4)**2+1.)**0.5)
	    ars2 = alog(x(j)*1.e+4+((x(j)*1.e+4)**2+1.)**0.5)
	    a    = (y(j)*x(j)*1.e+4-y(j-1)*x(j-1)*1.e+4)/(ars2-ars1)
	    b    = (y(j)*x(j)*1.e+4*ars1-y(j-1)*x(j-1)*1.e+4*ars2) / (ars1-ars2)
      yi   = a/xi/1.e+4*ars + b/xi/1.e+4
      return
   else if (xi == x(j)) then
      yi = y(j)
      return
   else
      j=j+1
   end if
   if (j <= n) then
      goto 30
   else
      yi = y(n)
   end if

end subroutine hyper_interp
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine CNinverse(Gmax,nspr,gamma,R,CNinv)
   
  implicit none

  real*8,    intent(IN)      :: Gmax
  integer, intent(IN)      :: nspr
  real,    intent(IN)      :: gamma(nspr)
  real,    intent(IN)      :: R(nspr)

  real          :: CNinv(nspr-1), sumdummy
  integer       :: i

  CNinv     = 0.0
  sumdummy  = 0.0
  do i = 1,nspr-1			
    CNinv(i) = ((gamma(i+1)/2.0- gamma(i)/2.0)/(R(i+1)-R(i)))- 0.5/Gmax- sumdummy
    sumdummy = sumdummy+ CNinv(i)
  enddo
 
end subroutine Cninverse
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
subroutine nonline_curve_hyper(Tdomain,Enonli,Gmax,isoil,R,H,CNinv)
   
  use sdomain
   
  implicit none

  type(nonli),    intent(INOUT),  target  :: Enonli
  type(domain),   intent(IN),     target  :: Tdomain
  integer,        intent(IN)              :: isoil   
  real(kind=8),   intent(IN)              :: Gmax
   
  integer   :: i,j,n			
  real      :: x0,xu,dx
  real      :: gamma(Enonli%Nspr)
  real      :: G(Enonli%Nspr) 
  real      :: R(Enonli%Nspr),H(Enonli%Nspr), CNinv(Enonli%Nspr-1)		
  integer   :: nspr

  nspr  = Enonli%Nspr

  x0    = -6.0		
  xu    = log10(0.10)			
  dx    = (xu-x0)/(nspr- 1)  


  do i = 1,nspr             
    gamma(i) = 10.0**(x0 + dx*(i- 1))
    G(i)     = (1.0/(1.0+ abs(gamma(i)/Enonli%gammaref(isoil))))
    R(i)     = Gmax* G(i)* gamma(i)
    write(33,*)   gamma(i), R(i), Enonli%gammaref(isoil)                       
  enddo

  call CNinverse(Gmax,nspr,gamma,R,CNinv)
 
  return
 
 end subroutine nonline_curve_hyper
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !Subroutine to find the inverse of a square matrix
        !Author : Louisda16th a.k.a Ashwith J. Rego
        !Reference : Algorithm has been well explained in:
        !http://math.uww.edu/~mcfarlat/inverse.htm
        !http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html

SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
      IMPLICIT NONE
      !Declarations
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
        REAL*8, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
        REAL*8, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

        LOGICAL :: FLAG = .TRUE.
        INTEGER :: i, j, k, l
        REAL*8 :: m
        REAL*8, DIMENSION(n,2*n) :: augmatrix !augmented matrix

!Augment input matrix with an identity matrix
        DO i = 1, n
      DO j = 1, 2*n
         IF (j <= n ) THEN
            augmatrix(i,j) = matrix(i,j)
         ELSE IF ((i+n) == j) THEN
            augmatrix(i,j) = 1.0D0
         Else
            augmatrix(i,j) = 0.0D0
         ENDIF
      END DO
        END DO

!Reduce augmented matrix to upper traingular form
        DO k =1, n-1
      IF (augmatrix(k,k) == 0.0D0) THEN
         FLAG = .FALSE.
         DO i = k+1, n
            IF (augmatrix(i,k) /= 0.0D0) THEN
               DO j = 1,2*n
                  augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
               END DO
               FLAG = .TRUE.
               EXIT
            ENDIF
            IF (FLAG .EQV. .FALSE.) THEN
               PRINT*, "Matrix is non - invertible"
               inverse = 0.0D0
               errorflag = -1
               return
            ENDIF
         END DO
      ENDIF
      DO j = k+1, n
         m = augmatrix(j,k)/augmatrix(k,k)
         DO i = k, 2*n
            augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
         END DO
      END DO
        END DO

!Test for invertibility
        DO i = 1, n
      IF (augmatrix(i,i) == 0.0D0) THEN
         PRINT*, "Matrix is non - invertible"
         inverse = 0.0D0
         errorflag = -1
         return
      ENDIF
        END DO

!Make diagonal elements as 1
        DO i = 1 , n
      m = augmatrix(i,i)
                DO j = i , (2 * n)
         augmatrix(i,j) = (augmatrix(i,j) / m)
                END DO
        END DO

!Reduced right side half of augmented matrix to identity matrix
        DO k = n-1, 1, -1
                DO i =1, k
         m = augmatrix(i,k+1)
         DO j = k, (2*n)
            augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
         END DO
      END DO
        END DO

!store answer
        DO i =1, n
      DO j = 1, n
         inverse(i,j) = augmatrix(i,j+n)
      END DO
        END DO
        errorflag = 0
END SUBROUTINE FINDinv
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module nonlinear
