!-------------------------Module de parametres globaux-----------------
!***********************************************************************
module gp
!***********************************************************************
  integer, parameter :: P8 = selected_real_kind(8)
  integer, parameter :: DP = KIND(1.0D0)
  integer, parameter :: SP = KIND(1.0)
  integer, parameter :: I4 = selected_int_kind(9)
!
  real(P8),parameter :: pi=3.141592653589793_P8
  real(P8),parameter :: pi2=2._DP*PI
!
!***********************************************************************
end module gp
!***********************************************************************
