! Module   : core_table
! Purpose  : monovariate tables

$include 'core.inc'

module core_table

  ! Uses

  use core_kinds
  use core_spline

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: deriv
  public :: integ
  public :: integ_def

  ! Procedures

contains

  function deriv (x, y, deriv_type, dy_dx_a, dy_dx_b) result (dy_dx)

    real(WP), intent(in)           :: x(:)
    real(WP), intent(in)           :: y(:)
    character(*), intent(in)       :: deriv_type
    real(WP), optional, intent(in) :: dy_dx_a
    real(WP), optional, intent(in) :: dy_dx_b
    real(WP)                       :: dy_dx(SIZE(x))

    type(spline_t) :: sp

    ! Calculate dy/dx via spline fitting

    sp = spline_t(x, y, deriv_type, dy_dx_a, dy_dx_b)

    dy_dx = sp%deriv()

    ! Finish

    return

  end function deriv

!****

  function integ (x, y, deriv_type, dy_dx_a, dy_dx_b) result (Y_)

    real(WP), intent(in)           :: x(:)
    real(WP), intent(in)           :: y(:)
    character(*), intent(in)       :: deriv_type
    real(WP), optional, intent(in) :: dy_dx_a
    real(WP), optional, intent(in) :: dy_dx_b
    real(WP)                       :: Y_(SIZE(x))
    
    type(spline_t) :: sp

    ! Calculate int_0^y y dx via spline fitting

    sp = spline_t(x, y, deriv_type, dy_dx_a, dy_dx_b)

    Y_ = sp%integ()

    ! Finish

    return

  end function integ

!****

  function integ_def (x, y, deriv_type, dy_dx_a, dy_dx_b) result (Y_def)

    real(WP), intent(in)           :: x(:)
    real(WP), intent(in)           :: y(:)
    character(*), intent(in)       :: deriv_type
    real(WP), optional, intent(in) :: dy_dx_a
    real(WP), optional, intent(in) :: dy_dx_b
    real(WP)                       :: Y_def
    
    real(WP) :: Y_(SIZE(x))

    ! Calculate int_{x_min}^{x_max} y dx via spline fitting

    Y_ = integ(x, y, deriv_type, dy_dx_a, dy_dx_b)

    Y_def = Y_(SIZE(x))

    ! Finish

    return

  end function integ_def

end module core_table
