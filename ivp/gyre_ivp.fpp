! Module   : gyre_ivp
! Purpose  : solve initial-value problems across single intervals
!
! Copyright 2013 Rich Townsend
!
! This file is part of GYRE. GYRE is free software: you can
! redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, version 3.
!
! GYRE is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

$include 'core.inc'

module gyre_ivp

  ! Uses

  use core_kinds

  use gyre_jacobian
  use gyre_ext_arith
  use gyre_ivp_magnus
  use gyre_ivp_findiff

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: solve
  public :: recon
  public :: abscissa

  ! Procedures

contains

  subroutine solve (solver_type, jc, omega, x_a, x_b, E_l, E_r, S, use_real)

    character(LEN=*), intent(in)     :: solver_type
    class(jacobian_t), intent(in)    :: jc
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_a
    real(WP), intent(in)             :: x_b
    complex(WP), intent(out)         :: E_l(:,:)
    complex(WP), intent(out)         :: E_r(:,:)
    type(ext_complex_t), intent(out) :: S
    logical, intent(in), optional    :: use_real

    ! Solve the IVP across the interval x_a -> x_b

    select case(solver_type)
    case('FINDIFF','FINDIFF_GL2')
       call solve_findiff_GL2(jc, omega, x_a, x_b, E_l, E_r, S)
    case('FINDIFF_GL4')
       call solve_findiff_GL4(jc, omega, x_a, x_b, E_l, E_r, S)
    case('MAGNUS_GL2')
       call solve_magnus_GL2(jc, omega, x_a, x_b, E_l, E_r, S, use_real)
    case('MAGNUS_GL4')
       call solve_magnus_GL4(jc, omega, x_a, x_b, E_l, E_r, S, use_real)
    case('MAGNUS_GL6')
       call solve_magnus_GL6(jc, omega, x_a, x_b, E_l, E_r, S, use_real)
    case default
       $ABORT(Invalid solver_type)
    end select

    ! Finish

    return

  end subroutine solve

!****

  subroutine recon (solver_type, jc, omega, x_a, x_b, y_a, y_b, x, y, use_real)

    character(LEN=*), intent(in)  :: solver_type
    class(jacobian_t), intent(in) :: jc
    complex(WP), intent(in)       :: omega
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    complex(WP), intent(in)       :: y_a(:)
    complex(WP), intent(in)       :: y_b(:)
    real(WP), intent(in)          :: x(:)
    complex(WP), intent(out)      :: y(:,:)
    logical, intent(in), optional :: use_real

    ! Reconstruct the IVP solution within the interval x_a -> x_b

    select case(solver_type)
    case('FINDIFF','FINDIFF_GL2')
       call recon_findiff_GL2(jc, omega, x_a, x_b, y_a, y_b, x, y)
    case('FINDIFF_GL4')
       call recon_findiff_GL4(jc, omega, x_a, x_b, y_a, y_b, x, y)
    case('MAGNUS_GL2')
       call recon_magnus_GL2(jc, omega, x_a, x_b, y_a, y_b, x, y, use_real)
    case('MAGNUS_GL4')
       call recon_magnus_GL4(jc, omega, x_a, x_b, y_a, y_b, x, y, use_real)
    case('MAGNUS_GL6')
       call recon_magnus_Gl6(jc, omega, x_a, x_b, y_a, y_b, x, y, use_real)
    case default
       $ABORT(Invalid solver_type)
    end select

    ! Finish

    return

  end subroutine recon

!****

  function abscissa (solver_type, x_a, x_b) result (x)

    character(LEN=*), intent(in) :: solver_type
    real(WP), intent(in)         :: x_a
    real(WP), intent(in)         :: x_b
    real(WP), allocatable        :: x(:)

    ! Determine the abscissa used for solving/reconstructing across
    ! the interval x_a -> x_b

    select case(solver_type)
    case('FINDIFF','FINDIFF_GL2')
       x = abscissa_findiff_GL2(x_a, x_b)
    case('FINDIFF_GL4')
       x = abscissa_findiff_GL4(x_a, x_b)
    case('MAGNUS_GL2')
       x = abscissa_magnus_GL2(x_a, x_b)
    case('MAGNUS_GL4')
       x = abscissa_magnus_GL4(x_a, x_b)
    case('MAGNUS_GL6')
       x = abscissa_magnus_GL6(x_a, x_b)
    case default
       $ABORT(Invalid solver_type)
    end select

    ! Finish

    return

  end function abscissa

end module gyre_ivp
