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

!****

  function deriv (x, y) result (dy_dx)

    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: y(:)
    real(WP)             :: dy_dx(SIZE(x))

    type(spline_t) :: sp

    ! Calculate dy/dx via spline fitting

    call sp%init(x, y)

    dy_dx = sp%deriv()

    ! Finish

    return

  end function deriv

!****

  function integ (x, y) result (Y_)

    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: y(:)
    real(WP)             :: Y_(SIZE(x))
    
    type(spline_t) :: sp

    ! Calculate int_0^y y dx via spline fitting

    call sp%init(x, y)

    Y_ = sp%integ()

    ! Finish

    return

  end function integ

!****

  function integ_def (x, y) result (Y_def)

    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: y(:)
    real(WP)             :: Y_def
    
    type(spline_t) :: sp

    ! Calculate int_{x_min}^{x_max} y dx via spline fitting

    call sp%init(x, y)

    Y_def = sp%integ_def()

    ! Finish

    return

  end function integ_def

end module core_table
