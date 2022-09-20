! Module   : gyre_parfait_model
! Purpose  : stellar parf (piecewise analytic representation) model
!
! Copyright 2022 Rich Townsend & The GYRE Team
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

module gyre_parfait_model

  ! Uses

  use core_kinds, only: WP

  use gyre_constants
  use gyre_grid
  use gyre_math
  use gyre_model
  use gyre_point
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (model_t) :: parfait_model_t
     private
     type(grid_t)          :: gr
     real(WP), allocatable :: m(:)
     real(WP), allocatable :: x(:)
     real(WP), allocatable :: y(:)
     real(WP), allocatable :: z(:)
     real(WP), allocatable :: alpha(:)
     real(WP), allocatable :: Gamma_1(:)
     real(WP)              :: Omega_rot
     integer               :: s_i
     integer               :: s_o
   contains
     private
     procedure, public :: coeff
     procedure         :: coeff_V_2_
     procedure         :: coeff_As_
     procedure         :: coeff_U_
     procedure         :: coeff_c_1_
     procedure         :: coeff_Gamma_1_
     procedure, public :: dcoeff
     procedure         :: dcoeff_V_2_
     procedure         :: dcoeff_As_
     procedure         :: dcoeff_U_
     procedure         :: dcoeff_c_1_
     procedure, public :: is_defined
     procedure, public :: is_vacuum
     procedure, public :: Delta_p
     procedure, public :: Delta_g
     procedure, public :: grid
  end type parfait_model_t

  ! Interfaces

  interface parfait_model_t
     module procedure parfait_model_t_
  end interface parfait_model_t

  ! Access specifiers

  private

  public :: parfait_model_t

  ! Procedures

contains

  function parfait_model_t_ (x, m, Gamma_1) result (ml)

    real(WP), intent(in)  :: x(:)
    real(WP), intent(in)  :: m(:)
    real(WP), intent(in)  :: Gamma_1(:)
    type(parfait_model_t) :: ml

    integer               :: n
    real(WP), allocatable :: y(:)
    real(WP), allocatable :: z(:)
    real(WP), allocatable :: alpha(:)
    integer               :: k
    real(WP)              :: beta
    real(WP), allocatable :: gr_x(:)

    $CHECK_BOUNDS(SIZE(m), SIZE(x))
    $CHECK_BOUNDS(SIZE(Gamma_1), SIZE(x)-1)

    ! Construct the parfait_model_t

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Constructing parfait (piecewise analytic representation) model'
100    format(A)
    endif

    ! Sanity checks

    n = SIZE(x)

    $ASSERT(x(1) == 0._WP,Radius coordinate must begin at center)
    $ASSERT(m(1) == 0._WP,Mass coordinate must begin at center)

    $ASSERT(ALL(x(2:n) > x(:n-1)),Non-monotonic radius coordinate)
    $ASSERT(ALL(m(2:n) > m(:n-1)),Non-monotonic mass coordinate)

    ! Integrate the dimensionless equation of hydrostatic equibrium,
    ! downward from the surface

    allocate(y(n))

    y(n) = 0._WP

    allocate(z(n-1))
    allocate(alpha(n-1))

    y_loop : do k = n-1, 1, -1

       ! Evaluate the dimensionless density z

       z(k) = 3._WP*(m(k+1) - m(k))/(4._WP*PI*(x(k+1)**3 - x(k)**3))

       ! Evaluate the dimensionless pressure y. alpha is a helper
       ! array used later to evaluate y between grid points

       alpha(k) = z(k)*(m(k+1) - m(k))*(x(k+1)**2 - x(k)**2)/(2._WP*(x(k+1)**3 - x(k)**3))

       if (k > 1) then

          beta = -z(k)*(m(k+1)*x(k)**3 - m(k)*x(k+1)**3)*(x(k)**(-1) - x(k+1)**(-1))/(x(k+1)**3 - x(k)**3)

          y(k) = y(k+1) + alpha(k) + beta

       else

          y(k) = y(k+1) + alpha(k)

       end if

    end do y_loop

    ! Store data

    ml%m = m
    ml%x = x
    ml%y = y
    ml%z = z
    ml%alpha = alpha

    ml%Gamma_1 = Gamma_1

    ! Set up the grid

    allocate(gr_x(2*n-2))

    gr_x(1) = x(1)
    gr_x(2:2*n-4:2) = x(2:n-1)
    gr_x(3:2*n-3:2) = x(2:n-1)
    gr_x(2*n-2) = x(n)

    ml%gr = grid_t(gr_x)

    ! Other stuff

    ml%Omega_rot = 0._WP

    ml%s_i = ml%gr%s_i()
    ml%s_o = ml%gr%s_o()

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Created', n, 'points'
110    format(3X,A,1X,I0,1X,A)
    endif

    ! Finish

    return

  end function parfait_model_t_

  !****

  function coeff (this, i, pt)

    class(parfait_model_t), intent(in) :: this
    integer, intent(in)                :: i
    type(point_t), intent(in)          :: pt
    real(WP)                           :: coeff

    $ASSERT_DEBUG(i >= 1 .AND. i <= I_LAST,Invalid index)
    $ASSERT_DEBUG(this%is_defined(i),Undefined coefficient)

    $ASSERT_DEBUG(pt%s >= this%s_i .AND. pt%s <= this%s_o,Invalid segment)

    ! Evaluate the i'th coefficient

    select case (i)
    case (I_V_2)
       coeff = this%coeff_V_2_(pt)
    case (I_AS)
       coeff = this%coeff_As_(pt)
    case (I_U)
       coeff = this%coeff_U_(pt)
    case (I_C_1)
       coeff = this%coeff_c_1_(pt)
    case (I_GAMMA_1)
       coeff = this%coeff_Gamma_1_(pt)
    case (I_OMEGA_ROT)
       coeff = this%Omega_rot
    end select

    ! Finish

    return

  end function coeff

  !****

  function coeff_V_2_ (this, pt) result (coeff)

    class(parfait_model_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    real(WP)                           :: coeff

    real(WP) :: w_m1
    real(WP) :: w_2
    real(WP) :: y

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),V_2 evaluation at vacuum point)

    ! Evaluate the V_2 coefficient

    associate (k => pt%s)

      ! Set up weight functions

      if (k > 1) then

         w_m1 = (pt%x**(-1) - this%x(k)**(-1))/(this%x(k+1)**(-1) - this%x(k)**(-1))

      else

         if (pt%x > 0._WP) then
            w_m1 = 1._WP
         else
            w_m1 = 0._WP
         endif

      endif
         
      w_2 = (pt%x**2 - this%x(k)**2)/(this%x(k+1)**2 - this%x(k)**2)

      ! Evaluate the dimensionless pressure
 
      y = (1._WP - w_m1)*this%y(k) + w_m1*this%y(k+1) - (w_2 - w_m1)*this%alpha(k)

      ! Evaluate the coefficient

      coeff = this%z(k)/(this%coeff_c_1_(pt)*y)

    end associate

    ! Finish

    return

  end function coeff_V_2_

  !****

  function coeff_As_ (this, pt) result (coeff)

    class(parfait_model_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    real(WP)                           :: coeff

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),As evaluation at vacuum point)

    ! Evaluate the As coefficient

    coeff = -this%coeff_V_2_(pt)*pt%x**2/this%coeff_Gamma_1_(pt)

    ! Finish

    return

  end function coeff_As_

  !****

  function coeff_U_ (this, pt) result (coeff)

    class(parfait_model_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    real(WP)                           :: coeff

    ! Evaluate the U coefficient

    associate (k => pt%s)

      coeff = 4._WP*PI*this%coeff_c_1_(pt)*this%z(k)

    end associate

    ! Finish

    return

  end function coeff_U_

  !****

  function coeff_c_1_ (this, pt) result (coeff)

    class(parfait_model_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    real(WP)                           :: coeff

    real(WP) :: w_3
    real(WP) :: m

    ! Evaluate the c_1 coefficient

    associate (k => pt%s)

      if (k > 1) then

         ! Set up weight function

         w_3 = (pt%x**3 - this%x(k)**3)/(this%x(k+1)**3 - this%x(k)**3)

         ! Evaluate the dimensionless mass

         m = (1._WP - w_3)*this%m(k) + w_3*this%m(k+1)

         ! Evaluate the coefficient

         coeff = pt%x**3/m

      else

         coeff = 3._WP/(4._WP*PI*this%z(1))

      end if

    end associate

    ! Finish

    return

  end function coeff_c_1_

  !****

  function coeff_Gamma_1_ (this, pt) result (coeff)

    class(parfait_model_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    real(WP)                           :: coeff

    ! Evaluate the Gamma_1 coefficient

    associate (k => pt%s)

      coeff = this%Gamma_1(k)

    end associate
    ! Finish

    return

  end function coeff_Gamma_1_

  !****

  function dcoeff (this, i, pt)

    class(parfait_model_t), intent(in) :: this
    integer, intent(in)                :: i
    type(point_t), intent(in)          :: pt
    real(WP)                           :: dcoeff

    $ASSERT_DEBUG(i >= 1 .AND. i <= I_LAST,Invalid index)
    $ASSERT_DEBUG(this%is_defined(i),Undefined coefficient)

    $ASSERT_DEBUG(pt%s >= this%s_i .AND. pt%s <= this%s_o,Invalid segment)

    ! Evaluate the i'th coefficient

    select case (i)
    case (I_V_2)
       dcoeff = this%dcoeff_V_2_(pt)
    case (I_AS)
       dcoeff = this%dcoeff_As_(pt)
    case (I_U)
       dcoeff = this%dcoeff_U_(pt)
    case (I_C_1)
       dcoeff = this%dcoeff_c_1_(pt)
    case (I_GAMMA_1)
       dcoeff = 0._WP
    case (I_OMEGA_ROT)
       dcoeff = 0._WP
    end select

    ! Finish

    return

  end function dcoeff

  !****

  function dcoeff_V_2_ (this, pt) result (dcoeff)

    class(parfait_model_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    real(WP)                           :: dcoeff

    ! Evaluate the logarithmic derivative of the V_2 coefficient

    dcoeff = this%coeff_V_2_(pt)*pt%x**2 + this%coeff_U_(pt) - 3._WP

    ! Finish

    return

  end function dcoeff_V_2_

  !****

  function dcoeff_As_ (this, pt) result (dcoeff)

    class(parfait_model_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    real(WP)                           :: dcoeff

    ! Evaluate the logarithmic derivative of the As coefficient

    dcoeff = -(this%dcoeff_V_2_(pt)*pt%x**2 + 2._WP*this%coeff_V_2_(pt)*pt%x)/this%coeff_Gamma_1_(pt)

    ! Finish

    return

  end function dcoeff_As_

  !****

  function dcoeff_U_ (this, pt) result (dcoeff)

    class(parfait_model_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    real(WP)                           :: dcoeff

    ! Evaluate the logarithmic derivative of the U coefficient, using
    ! eqn. (21) of Takata (2006) with dlnrho/dlnr = 0

    dcoeff = 3._WP - this%coeff_U_(pt)

    ! Finish

    return

  end function dcoeff_U_

  !****

  function dcoeff_c_1_ (this, pt) result (dcoeff)

    class(parfait_model_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    real(WP)                           :: dcoeff

    ! Evaluate the logarithmic derivative of the c_1 coefficient,
    ! using eqn. (20) of Takata (2006)

    dcoeff = 3._WP - this%coeff_U_(pt)

    ! Finish

    return

  end function dcoeff_c_1_

  !****

  function is_defined (this, i)

    class(parfait_model_t), intent(in) :: this
    integer, intent(in)                :: i
    logical                            :: is_defined

    $ASSERT_DEBUG(i >= 1 .AND. i <= I_LAST,Invalid index)

    ! Return the definition status of the i'th coefficient

    select case (i)
    case (I_V_2, I_AS, I_U, I_C_1, I_GAMMA_1, I_OMEGA_ROT)
       is_defined = .TRUE.
    case default
       is_defined = .FALSE.
    end select

    ! Finish

    return

  end function is_defined

  !****

  function is_vacuum (this, pt)

    class(parfait_model_t), intent(in) :: this
    type(point_t), intent(in)          :: pt
    logical                            :: is_vacuum

    $ASSERT_DEBUG(pt%s >= this%s_i .AND. pt%s <= this%s_o,Invalid segment)

    ! Return whether the point is a vacuum

    is_vacuum = (1._WP - pt%x**2) == 0._WP

    ! Finish

    return

  end function is_vacuum

  !****

  function Delta_p (this, x_i, x_o)

    class(parfait_model_t), intent(in) :: this
    real(WP), intent(in)               :: x_i
    real(WP), intent(in)               :: x_o
    real(WP)                           :: Delta_p

    ! Evaluate the dimensionless p-mode frequency separation

    $ABORT(Not yet implemented)

    ! Finish

    return

  end function Delta_p

  !****

  function Delta_g (this, x_i, x_o, lambda)

    class(parfait_model_t), intent(in) :: this
    real(WP), intent(in)               :: x_i
    real(WP), intent(in)               :: x_o
    real(WP), intent(in)               :: lambda
    real(WP)                           :: Delta_g

    ! Evaluate the dimensionless g-mode inverse period separation

    $ABORT(Not yet implemented)

    ! Finish

    return

  end function Delta_g

  !****

  function grid (this) result (gr)

    class(parfait_model_t), intent(in) :: this
    type(grid_t)                      :: gr

    ! Return the grid

    gr = this%gr

    ! Finish

    return

  end function grid

end module gyre_parfait_model
