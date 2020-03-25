! Module   : core_kinds
! Purpose  : kind types

module core_kinds

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: I4 = SELECTED_INT_KIND(9)
  integer, parameter :: I8 = SELECTED_INT_KIND(14)

  integer, parameter :: SP = KIND(0.)
  integer, parameter :: DP = KIND(0.D0)
  integer, parameter :: QP = SELECTED_REAL_KIND(32)

  $if($DOUBLE_PRECISION)
  integer, parameter :: WP = DP
  $else
  integer, parameter :: WP = SP
  $endif

  ! Acess specifiers

  private

  public :: I4
  public :: I8
  public :: SP
  public :: DP
  public :: WP
  public :: QP

end module core_kinds
