! Module   : gyre_oni_model
! Purpose  : stellar onion (piecewise constant density) model
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

module gyre_oni_model

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

  type, extends (model_t) :: oni_model_t
     private
     type(grid_t)          :: gr
     real(WP), allocatable :: m(:)
     real(WP), allocatable :: x(:)
     real(WP), allocatable :: y(:)
     real(WP), allocatable :: z(:)
     real(WP), allocatable :: alpha(:)
     real(WP), public      :: M_star
     real(WP), public      :: R_star
     real(WP)              :: Gamma_1
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
     procedure, public :: dcoeff
     procedure         :: dcoeff_V_2_
     procedure         :: dcoeff_As_
     procedure         :: dcoeff_U_
     procedure         :: dcoeff_c_1_
     procedure, public :: M_r
     procedure, public :: P
     procedure, public :: rho
     procedure, public :: is_defined
     procedure, public :: is_vacuum
     procedure, public :: Delta_p
     procedure, public :: Delta_g
     procedure, public :: grid
  end type oni_model_t

  ! Interfaces

  interface oni_model_t
     module procedure oni_model_t_
  end interface oni_model_t

  ! Access specifiers

  private

  public :: oni_model_t

  ! Procedures

contains

  function oni_model_t_ (M_r, r, Gamma_1) result (ml)

    real(WP), intent(in) :: M_r(:)
    real(WP), intent(in) :: r(:)
    real(WP), intent(in) :: Gamma_1
    type(oni_model_t)    :: ml

    integer               :: n
    real(WP), allocatable :: m(:)
    real(WP), allocatable :: x(:)
    real(WP), allocatable :: y(:)
    real(WP), allocatable :: z(:)
    real(WP), allocatable :: alpha(:)
    integer               :: k
    real(WP)              :: beta
    real(WP), allocatable :: gr_x(:)

    $CHECK_BOUNDS(SIZE(M_r), SIZE(r))

    ! Construct the oni_model_t

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Constructing onion (piecewise constant density) model'
100    format(A)
    endif

    ! Sanity checks

    n = SIZE(r)

    $ASSERT(ALL(M_r(2:) >= M_r(:n-1)),Non-monotonic mass coordinate)
    $ASSERT(ALL(r(2:) > r(:n-1)),Non-monotonic radial coordinate)

    $ASSERT(M_r(1) == 0._WP,Missing center)
    $ASSERT(r(1) == 0._WP,Missing center)

    ! Set stellar parameters

    ml%M_star = M_r(n)
    ml%R_star = r(n)

    ! Set up the dimensionless mass and radius coordinates

    m = M_r/ml%M_star
    x = r/ml%R_star

    ! Integrate the dimensionless equation of hydrostatic equibrium,
    ! downward from the surface

    allocate(y(n))

    y(n) = 0._WP

    allocate(z(n-1))
    allocate(alpha(n-1))

    y_loop : do k = n-1, 1, -1

       ! Evaluate the dimensionless density

       z(k) = 3._WP*(m(k+1) - m(k))/(4._WP*PI*(x(k+1)**3 - x(k)**3))

       ! Evaluate the dimensionless pressure

       alpha(k) = z(k)*(m(k+1) - m(k))*(x(k+1)**2 - x(k)**2)/(2._WP*(x(k+1)**3 - x(k)**3))

       if (x(k) /= 0._WP) then

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

    ! Set up the grid

    allocate(gr_x(2*n-2))

    gr_x(1) = x(1)
    gr_x(2:2*n-4:2) = x(2:n-1)
    gr_x(3:2*n-3:2) = x(2:n-1)
    gr_x(2*n-2) = x(n)

    ml%gr = grid_t(gr_x)

    ! Other stuff

    ml%Gamma_1 = Gamma_1

    ml%Omega_rot = 0._WP

    ml%s_i = ml%gr%s_i()
    ml%s_o = ml%gr%s_o()

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Created', n, 'points'
110    format(3X,A,1X,I0,1X,A)
    endif

    ! Finish

    return

  end function oni_model_t_

  !****

  function coeff (this, i, pt)

    class(oni_model_t), intent(in) :: this
    integer, intent(in)            :: i
    type(point_t), intent(in)      :: pt
    real(WP)                       :: coeff

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
       coeff = this%Gamma_1
    case (I_OMEGA_ROT)
       coeff = this%Omega_rot
    end select

    ! Finish

    return

  end function coeff

  !****

  function coeff_V_2_ (this, pt) result (coeff)

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: coeff

    real(WP) :: w_m1
    real(WP) :: w_2
    real(WP) :: y

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),V_2 evaluation at vacuum point)

    ! Evaluate the V_2 coefficient

    if (pt%x == 0._WP) then

       coeff = this%z(1)/(this%coeff_c_1_(pt)*this%y(1))

    else

       associate(k => pt%s)

         ! Set up weight functions

         if (this%x(k) == 0._WP) then
            w_m1 = 1._WP
         else
            w_m1 = (pt%x**(-1) - this%x(k)**(-1))/(this%x(k+1)**(-1) - this%x(k)**(-1))
         endif
         
         w_2 = (pt%x**2 - this%x(k)**2)/(this%x(k+1)**2 - this%x(k)**2)

         ! Evaluate the dimensionless pressure
 
         y = (1._WP - w_m1)*this%y(k) + w_m1*this%y(k+1) - (w_2 - w_m1)*this%alpha(k)

         ! Evaluate the coefficient

         coeff = this%z(k)/(this%coeff_c_1_(pt)*y)

       end associate

    end if
         
    ! Finish

    return

  end function coeff_V_2_

  !****

  function coeff_As_ (this, pt) result (coeff)

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: coeff

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),As evaluation at vacuum point)

    ! Evaluate the As coefficient

    coeff = -this%coeff_V_2_(pt)*pt%x**2/this%Gamma_1

    ! Finish

    return

  end function coeff_As_

  !****

  function coeff_U_ (this, pt) result (coeff)

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: coeff

    ! Evaluate the U coefficient

    associate(k => pt%s)

      coeff = 4._WP*PI*this%coeff_c_1_(pt)*this%z(k)

    end associate

    ! Finish

    return

  end function coeff_U_

  !****

  function coeff_c_1_ (this, pt) result (coeff)

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: coeff

    real(WP) :: w_3
    real(WP) :: m

    ! Evaluate the c_1 coefficient

    if (pt%x /= 0._WP) then

       associate(k => pt%s)

         ! Set up weight functions

         w_3 = (pt%x**3 - this%x(k)**3)/(this%x(k+1)**3 - this%x(k)**3)

         ! Evaluate the dimensionless mass

         m = (1._WP - w_3)*this%m(k) + w_3*this%m(k+1)

         ! Evaluate the coefficient

         coeff = pt%x**3/m

       end associate

    else

       coeff = this%x(2)**3/this%m(2)

    end if

    ! Finish

    return

  end function coeff_c_1_

  !****

  function dcoeff (this, i, pt)

    class(oni_model_t), intent(in) :: this
    integer, intent(in)            :: i
    type(point_t), intent(in)      :: pt
    real(WP)                       :: dcoeff

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

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: dcoeff

    ! Evaluate the logarithmic derivative of the V_2 coefficient

    dcoeff = this%coeff_V_2_(pt)*pt%x**2 + this%coeff_U_(pt) - 3._WP

    ! Finish

    return

  end function dcoeff_V_2_

  !****

  function dcoeff_As_ (this, pt) result (dcoeff)

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: dcoeff

    ! Evaluate the logarithmic derivative of the As coefficient

    dcoeff = -(this%dcoeff_V_2_(pt)*pt%x**2 + 2._WP*this%coeff_V_2_(pt)*pt%x)/this%Gamma_1

    ! Finish

    return

  end function dcoeff_As_

  !****

  function dcoeff_U_ (this, pt) result (dcoeff)

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: dcoeff

    ! Evaluate the logarithmic derivative of the U coefficient, using
    ! eqn. (21) of Takata (2006) with dlnrho/dlnr = 0

    dcoeff = 3._WP - this%coeff_U_(pt)

    ! Finish

    return

  end function dcoeff_U_

  !****

  function dcoeff_c_1_ (this, pt) result (dcoeff)

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: dcoeff

    ! Evaluate the logarithmic derivative of the c_1 coefficient,
    ! using eqn. (20) of Takata (2006)

    dcoeff = 3._WP - this%coeff_U_(pt)

    ! Finish

    return

  end function dcoeff_c_1_

  !****

  function M_r (this, pt)

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: M_r

    ! Evaluate the fractional mass coordinate

    M_r = this%M_star*(pt%x**3/this%coeff(I_C_1, pt))

    ! Finish

    return

  end function M_r
    
  !****

  function P (this, pt)

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: P

    ! Evaluate the total pressure

    P = (G_GRAVITY*this%M_star**2/(4._WP*PI*this%R_star**4))* &
        (this%coeff(I_U, pt)/(this%coeff(I_C_1, pt)**2*this%coeff(I_V_2, pt)))

    ! Finish

    return

  end function P
    
  !****

  function rho (this, pt)

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: rho

    ! Evaluate the density

    rho = (this%M_star/(4._WP*PI*this%R_star**3))*(this%coeff(I_U, pt)/this%coeff(I_C_1, pt))

    ! Finish

    return

  end function rho
    
  !****

  function is_defined (this, i)

    class(oni_model_t), intent(in) :: this
    integer, intent(in)            :: i
    logical                        :: is_defined

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

    class(oni_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    logical                        :: is_vacuum

    $ASSERT_DEBUG(pt%s >= this%s_i .AND. pt%s <= this%s_o,Invalid segment)

    ! Return whether the point is a vacuum

    is_vacuum = (1._WP - pt%x**2) == 0._WP

    ! Finish

    return

  end function is_vacuum

  !****

  function Delta_p (this, x_i, x_o)

    class(oni_model_t), intent(in) :: this
    real(WP), intent(in)           :: x_i
    real(WP), intent(in)           :: x_o
    real(WP)                       :: Delta_p

    ! Evaluate the dimensionless p-mode frequency separation

    Delta_p = 0.5_WP/(sqrt(2._WP/this%Gamma_1)*(asin(x_o)-asin(x_i)))

    ! Finish

    return

  end function Delta_p

  !****

  function Delta_g (this, x_i, x_o, lambda)

    class(oni_model_t), intent(in) :: this
    real(WP), intent(in)           :: x_i
    real(WP), intent(in)           :: x_o
    real(WP), intent(in)           :: lambda
    real(WP)                       :: Delta_g

    ! Evaluate the dimensionless g-mode inverse period separation

    Delta_g = 0._WP

    ! Finish

    return

  end function Delta_g

  !****

  function grid (this) result (gr)

    class(oni_model_t), intent(in) :: this
    type(grid_t)                   :: gr

    ! Return the grid

    gr = this%gr

    ! Finish

    return

  end function grid

end module gyre_oni_model
