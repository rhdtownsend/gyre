! Module   : test_gyre_func
! Purpose  : unit tests for gyre_func
!
! Copyright 2019 The GYRE Team
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

program test_gyre_func

  ! Uses

  use core_kinds
  use core_constants

  use gyre_func

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  logical :: status

  ! Run tests

  status = .TRUE.

  call test_factorial(status)
  call test_double_factorial(status)
  call test_legendre_P(status)
  call test_spherical_Y(status)

  ! Finish
  
contains

  subroutine test_factorial (status)

    logical, intent(inout) :: status

    real(WP), parameter :: TOL = 2._WP*EPSILON(1._WP)

    ! Run checks

    call check_result_r('factorial (n=0)', factorial(0), REAL(1, WP), 0._WP, status)
    call check_result_r('factorial (n=1)', factorial(1), REAL(1, WP), 0._WP, status)
    call check_result_r('factorial (n=2)', factorial(2), REAL(2, WP), 0._WP, status)
    call check_result_r('factorial (n=3)', factorial(3), REAL(6, WP), 0._WP, status)
    call check_result_r('factorial (n=4)', factorial(4), REAL(24, WP), 0._WP, status)
    call check_result_r('factorial (n=5)', factorial(5), REAL(120, WP), 0._WP, status)
    call check_result_r('factorial (n=6)', factorial(6), REAL(720, WP), 0._WP, status)
    call check_result_r('factorial (n=7)', factorial(7), REAL(5040, WP), 0._WP, status)
    call check_result_r('factorial (n=8)', factorial(8), REAL(40320, WP), 0._WP, status)
    call check_result_r('factorial (n=9)', factorial(9), REAL(362880, WP), 0._WP, status)
    call check_result_r('factorial (n=10)', factorial(10), REAL(3628800, WP), 0._WP, status)

  end subroutine test_factorial

  !****

  subroutine test_double_factorial (status)

    logical, intent(inout) :: status

    real(WP), parameter :: TOL = 0._WP*EPSILON(1._WP)

    ! Run checks

    call check_result_r('double factorial (n=0)', double_factorial(0), REAL(1, WP), 0._WP, status)
    call check_result_r('double factorial (n=1)', double_factorial(1), REAL(1, WP), 0._WP, status)
    call check_result_r('double factorial (n=2)', double_factorial(2), REAL(2, WP), 0._WP, status)
    call check_result_r('double factorial (n=3)', double_factorial(3), REAL(3, WP), 0._WP, status)
    call check_result_r('double factorial (n=4)', double_factorial(4), REAL(8, WP), 0._WP, status)
    call check_result_r('double factorial (n=5)', double_factorial(5), REAL(15, WP), 0._WP, status)
    call check_result_r('double factorial (n=6)', double_factorial(6), REAL(48, WP), 0._WP, status)
    call check_result_r('double factorial (n=7)', double_factorial(7), REAL(105, WP), 0._WP, status)
    call check_result_r('double factorial (n=8)', double_factorial(8), REAL(384, WP), 0._WP, status)
    call check_result_r('double factorial (n=9)', double_factorial(9), REAL(945, WP), 0._WP, status)
    call check_result_r('double factorial (n=10)', double_factorial(10), REAL(3840, WP), 0._WP, status)

  end subroutine test_double_factorial

  !****

  subroutine test_legendre_P (status)

    logical, intent(inout) :: status

    real(WP), parameter :: TOL = 4._WP*EPSILON(0._WP)

    integer, parameter :: N_X = 101
    real(WP), parameter :: X_MIN = -1._WP
    real(WP), parameter :: X_MAX = 1._WP

    integer  :: i
    real(WP) :: x

    do i = 1, N_X

       x = (X_MIN*(N_X - i) + X_MAX*(i - 1._WP))/(N_X - 1._WP)

       call check_result_r('legendre_P (0,0)', legendre_P(0, 0, x), 1._WP, TOL, status)

       call check_result_r('legendre_P (1,-1)', legendre_P(1, -1, x), SQRT(1._WP - x**2)/2._WP, TOL, status)
       call check_result_r('legendre_P (1,0)', legendre_P(1, 0, x), x, TOL, status)
       call check_result_r('legendre_P (1,1)', legendre_P(1, 1, x), -SQRT(1._WP - x**2), TOL, status)
       
       call check_result_r('legendre_P (2,-2)', legendre_P(2, -2, x), (1._WP - x**2)/8._WP, TOL, status)
       call check_result_r('legendre_P (2,-1)', legendre_P(2, -1, x), x*SQRT(1._WP - x**2)/2._WP, TOL, status)
       call check_result_r('legendre_P (2,0)', legendre_P(2, 0, x), (3._WP*x**2 - 1._WP)/2._WP, TOL, status)
       call check_result_r('legendre_P (2,1)', legendre_P(2 ,1, x), -3._WP*x*SQRT(1._WP - x**2), TOL, status)
       call check_result_r('legendre_P (2,2)', legendre_P(2 ,2 ,x), 3._WP*(1._WP - x**2), TOL, status)
       
    end do

  end subroutine test_legendre_P

  !****

  subroutine test_spherical_Y (status)

    logical, intent(inout) :: status

    real(WP), parameter    :: TOL = EPSILON(0._WP)
    real(WP), parameter    :: PHI = 0.0_WP
    complex(WP), parameter :: II = CMPLX(0._WP, 1._WP, KIND=WP)

    integer, parameter :: N_X = 1001
    real(WP), parameter :: X_MIN = -1._WP
    real(WP), parameter :: X_MAX = 1._WP

    integer  :: i
    real(WP) :: x
    real(WP) :: theta

    do i = 1, N_X

       x = (X_MIN*(N_X - i) + X_MAX*(i - 1._WP))/(N_X - 1._WP)

       theta = ACOS(x)

       call check_result_c('spherical_Y (0,0)', spherical_Y(0, 0, theta, PHI), &
            0.5_WP/SQRT(PI)*EXP(0*II*PHI), TOL, status)

       call check_result_c('spherical_Y (1,-1)', spherical_Y(1, -1, theta, PHI), &
            0.5_WP*SQRT(3._WP/(2._WP*PI))*SIN(theta)*EXP(-1*II*PHI), TOL, status)
       call check_result_c('spherical_Y (1,0)', spherical_Y(1, 0, theta, PHI), &
            0.5_WP*SQRT(3._WP/PI)*COS(theta)*EXP(0*II*PHI), TOL, status)
       call check_result_c('spherical_Y (1,1)', spherical_Y(1, 1, theta, PHI), &
            -0.5_WP*SQRT(3._WP/(2._WP*PI))*SIN(theta)*EXP(1*II*PHI), TOL, status)
       
       call check_result_c('spherical_Y (2,-2)', spherical_Y(2, -2, theta, PHI), &
            0.25_WP*SQRT(15._WP/(2._WP*PI))*SIN(theta)**2*EXP(-2*II*PHI), TOL, status)
       call check_result_c('spherical_Y (2,-1)', spherical_Y(2, -1, theta, PHI), &
            0.5_WP*SQRT(15._WP/(2._WP*PI))*SIN(theta)*COS(theta)*EXP(-1*II*PHI), TOL, status)
       call check_result_c('spherical_Y (2,0)', spherical_Y(2, 0, theta, PHI), &
            0.25_WP*SQRT(5._WP/PI)*(3._WP*COS(theta)**2 - 1._WP)*EXP(0*II*PHI), TOL, status)
       call check_result_c('spherical_Y (2,1)', spherical_Y(2, 1, theta, PHI), &
            -0.5_WP*SQRT(15._WP/(2._WP*PI))*SIN(theta)*COS(theta)*EXP(1*II*PHI), TOL, status)
       call check_result_c('spherical_Y (2,2)', spherical_Y(2, 2, theta, PHI), &
            0.25_WP*SQRT(15._WP/(2._WP*PI))*SIN(theta)**2*EXP(2*II*PHI), TOL, status)
       
    end do

  end subroutine test_spherical_Y

  !****

  subroutine check_result_r (label, a, b, tol, status)

    character(*), intent(in) :: label
    real(WP), intent(in)     :: a
    real(WP), intent(in)     :: b
    real(WP), intent(in)     :: tol
    logical, intent(inout)   :: status

    if (ABS(a - b) > tol) then
       write(*, 100) 'Failed in ', label, ', difference = ', a, b, ulp_diff(a, b)
100    format(A, A, A, E24.17, 1X, E24.17, 1X, I0)
       status = .FALSE.
    end if

    ! Finish

    return

  end subroutine check_result_r

  !****

  subroutine check_result_c (label, a, b, tol, status)

    character(*), intent(in) :: label
    complex(WP), intent(in)  :: a
    complex(WP), intent(in)  :: b
    real(WP), intent(in)     :: tol
    logical, intent(inout)   :: status

    if (ABS(REAL(a - b)) > tol) then
       write(*, 100) 'Failed in ', label, ', re difference = ', REAL(a), REAL(b)
100    format(A, A, A, E24.17, 1X, E24.17, 1X, I0)
       status = .FALSE.
    end if

    if (ABS(AIMAG(a - b)) > tol) then
       write(*, 100) 'Failed in ', label, ', im difference = ', AIMAG(a), AIMAG(b)
       status = .FALSE.
    end if

    ! Finish

    return

  end subroutine check_result_c

  !****

  function ulp_diff (a, b) result (u)

    real(WP), intent(in) :: a
    real(WP), intent(in) :: b
    integer(I8)          :: u

    real(WP) :: s
    real(WP) :: c

    ! Calculate by how many ulps b differs from a

    s = SIGN(1._WP, b-a)

    u = 0
    c = a

    search_loop : do
       if (c == b) exit search_loop
       c = NEAREST(c, s)
       u = u + 1
    end do search_loop

    ! Finish

    return

  end function ulp_diff

end program test_gyre_func
