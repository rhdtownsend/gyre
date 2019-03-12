! Module   : gyre_twopt_model
! Purpose  : stellar two-point model
!
! Copyright 2018 Rich Townsend
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

module gyre_twopt_model

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_grid
  use gyre_model
  use gyre_model_par
  use gyre_model_util
  use gyre_point

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (model_t) :: twopt_model_t
     private
     type(grid_t) :: gr
     real(WP)     :: V_o
     real(WP)     :: dV_o
     real(WP)     :: U_o
     real(WP)     :: dU_o
     real(WP)     :: c_1_i
     real(WP)     :: Gamma_1
     real(WP)     :: Omega_rot
     integer      :: s_i
     integer      :: s_o
   contains
     private
     procedure, public :: coeff
     procedure         :: coeff_V_2_
     procedure         :: coeff_As_
     procedure         :: coeff_U_
     procedure         :: coeff_c_1_
     procedure, public :: dcoeff
     procedure         :: dcoeff_V_2_
     procedure         :: rho_
     procedure         :: m_r_
     procedure         :: P_
     procedure         :: dcoeff_U_
     procedure         :: dcoeff_c_1_
     procedure, public :: is_defined
     procedure, public :: is_vacuum
     procedure, public :: Delta_p
     procedure, public :: Delta_g
     procedure, public :: grid
  end type twopt_model_t

  ! Interfaces

  interface twopt_model_t
     module procedure twopt_model_t_
  end interface twopt_model_t

  
  ! Access specifiers

  private

  public :: twopt_model_t

  ! Procedures

contains

  function twopt_model_t_ (ml_p, V_o, dV_o, U_o, dU_o, c_1_i) result (ml)

    use gyre_grid_weights
  
    type(model_par_t), intent(in) :: ml_p
    real(WP), intent(in)          :: V_o
    real(WP), intent(in)          :: dV_o
    real(WP), intent(in)          :: U_o
    real(WP), intent(in)          :: dU_o
    real(WP), intent(in)          :: c_1_i
    type(twopt_model_t)           :: ml

    real(WP), allocatable :: w(:)

    ! Construct the twopt_model_t

    ml%gr = grid_t([0._WP, 1._WP])

    ml%V_o = V_o
    ml%dV_o = dV_o
 
    ml%U_o = U_o
    ml%dU_o = dU_o

    ml%c_1_i = c_1_i

    ml%Gamma_1 = ml_p%Gamma_1

    if (ml_p%uniform_rot) then
       ml%Omega_rot = uniform_Omega_rot(ml_p)
    else
       ml%Omega_rot = 0._WP
    endif

    ml%s_i = ml%gr%s_i()
    ml%s_o = ml%gr%s_o()

    ! Finish

    return

  end function twopt_model_t_

  !****

  function coeff (this, i, pt)

    class(twopt_model_t), intent(in) :: this
    integer, intent(in)              :: i
    type(point_t), intent(in)        :: pt
    real(WP)                         :: coeff

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
    case (I_DELTA)
       coeff = 1._WP
    case (I_NABLA_AD)
       coeff = 0.4_WP
    case (I_OMEGA_ROT)
       coeff = this%Omega_rot
    end select

    ! Finish

    return

  end function coeff

  !****

  function coeff_V_2_ (this, pt) result (coeff)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: coeff
    real(WP)                         :: n_n
    real(WP)                         :: m_n
    real(WP)                         :: alpha
    real(WP)                         :: beta

    ! Evaluate the V_2 coefficient

    associate ( &
         x => pt%x, &
         m_n => 0._WP, &
         n_n => 0._WP, &
         alpha => 1.15927_WP, &
         beta => 317.76_WP)
   
    IF (is_vacuum(this,pt) .eqv. .FALSE.) THEN
      coeff = ((alpha * beta)/(60._WP*this%P_(pt)))*(x-1._WP)*(-1._WP - (1._WP + m_n)*x + (2._WP + m_n + n_n)&
        *x**2)*(20._WP + 15._WP*m_n*x - 12._WP*(3._WP + 2._WP*m_n + n_n)*x**2 + 10._WP * &
        (2._WP + m_n + n_n)*x**3)
    ELSE IF (is_vacuum(this,pt) .eqv. .TRUE.) THEN
        coeff = HUGE(0._WP)
    END IF
    end associate

    ! Finish

    return

  end function coeff_V_2_

  !****

  function coeff_As_ (this, pt) result (coeff)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: coeff

    ! Evaluate the As coefficient

    coeff = 3._WP - this%coeff_U_(pt) - this%coeff_V_2_(pt)*pt%x**2/this%Gamma_1 - this%dcoeff_U_(pt)

    ! Finish

    return

  end function coeff_As_

  !****

  function coeff_U_ (this, pt) result (coeff)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: coeff
    real(WP)                         :: n_n
    real(WP)                         :: m_n
    real(WP)                         :: alpha
    real(WP)                         :: beta

    ! Evaluate the U coefficient

    associate ( &
         x => pt%x, &
         m_n => 0._WP, &
         n_n => 0._WP, &
         alpha => 1.15927_WP, &
         beta => 317.76_WP)

      coeff = (60._WP*(1._WP + x*(m_n - (3._WP + 2*m_n + n_n)*x + (2._WP + m_n + n_n)*x**2)))&
        /(20._WP + 15._WP*m_n*x - 12._WP*(3._WP + 2._WP*m_n + n_n)*x**2 + 10._WP*(2._WP + m_n&
        + n_n)*x**3)

    end associate

    ! Finish

    return

  end function coeff_U_

  !****

  function coeff_c_1_ (this, pt) result (coeff)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: coeff

    ! Evaluate the c_1 coefficient

    associate ( &
         x => pt%x, &
         m_n => 0._WP, &
         n_n => 0._WP, &
         alpha => 1.15927_WP, &
         beta => 317.76_WP)
      
      !THIS NEED TOTAL MASS
      coeff = (4._WP + m_n + 2._WP*n_n)/(20._WP + 15._WP*m_n*x - 12._WP*(3._WP + 2._WP*m_n + n_n)*x**2 +&
        10._WP*(2._WP + m_n + n_n)*x**3)

    end associate

    ! Finish

    return

  end function coeff_c_1_

  !****

  function dcoeff (this, i, pt)

    class(twopt_model_t), intent(in) :: this
    integer, intent(in)              :: i
    type(point_t), intent(in)        :: pt
    real(WP)                         :: dcoeff

    $ASSERT_DEBUG(i >= 1 .AND. i <= I_LAST,Invalid index)
    $ASSERT_DEBUG(this%is_defined(i),Undefined coefficient)

    $ASSERT_DEBUG(pt%s >= this%s_i .AND. pt%s <= this%s_o,Invalid segment)

    ! Evaluate the i'th coefficient

    select case (i)
    case (I_V_2)
       dcoeff = this%dcoeff_V_2_(pt)
    case (I_AS)
       $ABORT(dcoeff for As not implemented)
    case (I_U)
       dcoeff = this%dcoeff_U_(pt)
    case (I_C_1)
       dcoeff = this%dcoeff_c_1_(pt)
    case (I_GAMMA_1)
       dcoeff = 0._WP
    case (I_DELTA)
       dcoeff = 0._WP
    case (I_NABLA_AD)
       dcoeff = 0._WP
    case (I_OMEGA_ROT)
       dcoeff = 0._WP
    end select

    ! Finish

    return

  end function dcoeff

  !****

  function dcoeff_V_2_ (this, pt) result (dcoeff)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: dcoeff
    real(WP)                         :: h_n
    real(WP)                         :: g_n
    real(WP)                         :: n_n
    real(WP)                         :: m_n
    real(WP)                         :: alpha
    real(WP)                         :: beta

    ! Evaluate the logarithmic derivative of the V_2 coefficient

    associate ( &
         x => pt%x, &
         m_n => 0._WP, &
         n_n => 0._WP, &
         alpha => 1.15927_WP, &
         beta => 317.76_WP)

      g_n = 2._WP*(336._WP*m_n + 34*m_n**2 - 68._WP*m_n*n_n + 45._WP*n_n**2 + 3._WP*(34._WP*m_n**2 &
        + m_n*(336._WP - 68._WP*n_n) + 15._WP*n_n*(-20._WP + 3._WP*n_n))*x + 6._WP*(m_n*(-889._WP &
        + 34._WP*m_n) - 68._WP*m_n*n_n + 45._WP*n_n**2 - 300._WP*(4._WP + n_n))*x**2 - &
        10._WP*(-396._WP + 281._WP*m_n**2 - 3._WP*n_n*(124._WP + 15._WP*n_n) + 4._WP*m_n*(70._WP&
        + 17._WP*n_n))*x**3 + 75._WP*(2._WP + m_n + n_n)*(62._WP + 53._WP*m_n + 9._WP*n_n)*x**4&
        - 1575._WP*(2._WP + m_n + n_n)**2*x**5 - 300._WP*(-3._WP + n_n + 5._WP*x))
      
      h_n = -900._WP - 336._WP*m_n - 34._WP*m_n**2 + 300._WP*n_n + 68._WP*m_n*n_n - 45._WP*n_n**2&
        - (900._WP + 34._WP*m_n**2 + m_n*(336._WP - 68._WP*n_n) + 15._WP*n_n*(-20._WP + 3._WP*n_n))*x &
        - 3._WP*(-500._WP + 336._WP*m_n + 34._WP*m_n**2 - 300._WP*n_n - 68._WP*m_n*n_n &
        + 45._WP*n_n**2)*x**2 + 4._WP*(-34._WP*m_n**2 - 45._WP*n_n**2 + 300._WP*(4._WP + n_n)&
        + m_n*(889._WP + 68._WP*n_n))*x**3 + 5._WP*(-396._WP + 281._WP*m_n**2 - 3._WP*n_n*(124._WP&
        + 15._WP*n_n) + 4._WP*m_n*(70._WP + 17._WP*n_n))*x**4 - 30._WP*(2._WP + m_n + n_n)*(62._WP&
        + 53._WP*m_n + 9._WP*n_n)*x**5 + 525._WP*(2._WP + m_n + n_n)**2*x**6

      dcoeff = x*((1._WP)/(x-1._WP) - (1._WP + m_n - 2._WP*(2._WP + m_n + n_n)*x)/(-1._WP &
        - (1._WP + m_n)*x + (2._WP + m_n + n_n)*x**2) + (15._WP*m_n - 24._WP*(3._WP + 2._WP*m_n &
        + n_n)*x + 30._WP*(2._WP + m_n + n_n)*x**2)/(20._WP + 15._WP*m_n*x - 12._WP*(3._WP + &
        2._WP*m_n + n_n)*x**2 + 10._WP*(2._WP + m_n + n_n)*x**3) + g_n/h_n)

    end associate

    ! Finish

    return

  end function dcoeff_V_2_

  !****

  function rho_ (this, pt) result (rho)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: rho
    real(WP)                         :: n_n
    real(WP)                         :: m_n
    real(WP)                         :: alpha
    real(WP)                         :: beta

    ! Evaluate the polynomial for density function

    associate ( &
         x => pt%x, &
         m_n => 0._WP, &
         n_n => 0._WP, &
         alpha => 1.15927_WP, &
         beta => 317.76_WP)

      rho = (2._WP + m_n + n_n)*x**3 + (-3._WP - 2._WP*m_n - n_n)*x**2 + m_n*x + 1._WP

    end associate

  end function rho_

  !****

  function m_r_ (this, pt) result (m_r)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: rho
    real(WP)                         :: m_r
    real(WP)                         :: m_n
    real(WP)                         :: n_n
    real(WP)                         :: alpha
    real(WP)                         :: beta

    ! Evaluate the polynomial for m_n profile

    associate ( &
         x => pt%x, &
         m_n => 0._WP, &
         n_n => 0._WP, &
         alpha => 1.15927_WP, &
         beta => 317.76_WP)
         
      m_r = x**3 * (beta/60._WP)*(20._WP + 15._WP*m_n*x - 12._WP*(3._WP + 2._WP*m_n&
        + n_n)*x**2 + 10._WP*(2._WP + m_n + n_n)*x**3)

    end associate

  end function m_r_

  !****

  function P_ (this, pt) result (Pressure)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: rho
    real(WP)                         :: Pressure
    real(WP)                         :: m_n
    real(WP)                         :: n_n
    real(WP)                         :: alpha
    real(WP)                         :: beta

    ! Evaluate the polynomial for Pressure profile  

    associate ( &
         x => pt%x, &
         m_n => 0._WP, &
         n_n => 0._WP, &
         alpha => 1.15927_WP, &
         beta => 317.76_WP)

      Pressure = -((alpha * beta)/(25200._WP))*(-336._WP*m_n - 34._WP*m_n**2 + 68._WP*m_n*n_n &
        - 45._WP*n_n**2 + 4900._WP*m_n*x**3 + 105._WP*(m_n*(-64._WP + 15._WP*m_n) - 32*(3._WP &
        + n_n))*x**4 - 252._WP*(-10._WP*(2._WP + n_n) + m_n*(17._WP + 18._WP*m_n + 9._WP*n_n))*x**5 &
        + 70._WP*(73._WP*m_n**2 + 12._WP*(3._WP + n_n)**2 + m_n*(194._WP + 73._WP*n_n))*x**6 - &
        1320._WP*(2._WP + m_n + n_n)*(3._WP + 2._WP*m_n + n_n)*x**7 + 525._WP*(2._WP + m_n + &
        n_n)**2*x**8 + 300._WP*(-3._WP + n_n + 14._WP*x**2))

    end associate

  end function P_

  !****

  function dcoeff_U_ (this, pt) result (dcoeff)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: dcoeff
    real(WP)                         :: m_n
    real(WP)                         :: n_n
    real(WP)                         :: alpha
    real(WP)                         :: beta

    ! Evaluate the logarithmic derivative of the U coefficient

    associate ( &
         x => pt%x, &
         m_n => 0._WP, &
         n_n => 0._WP, &
         alpha => 1.15927_WP, &
         beta => 317.76_WP)
    
      dcoeff = (1._WP/(x-1))+((2._WP + x + m_n*x)/((2._WP + m_n + n_n)*x**2 - (1._WP + m_n)*x - 1._WP) &
        + (60._WP + 6._WP*x*(5._WP*m_n - 2._WP*(3._WP + 2._WP*m_n + n_n)*x))/(20._WP + &
        15._WP*m_n*x - 12._WP*(3._WP + 2._WP*m_n + n_n)*x**2 + 10._WP*(2._WP + m_n + n_n)*x**3))

    end associate

    ! Finish

    return

  end function dcoeff_U_

  !****

  function dcoeff_c_1_ (this, pt) result (dcoeff)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: dcoeff

    ! Evaluate the logarithmic derivative of the c_1 coefficient

    dcoeff = 3._WP - this%coeff_U_(pt)

    ! Finish

    return

  end function dcoeff_c_1_

  !****

  function is_defined (this, i)

    class(twopt_model_t), intent(in) :: this
    integer, intent(in)              :: i
    logical                          :: is_defined

    $ASSERT_DEBUG(i >= 1 .AND. i <= I_LAST,Invalid index)

    ! Return the definition status of the i'th coefficient

    select case (i)
    case (I_V_2, I_AS, I_U, I_C_1, I_GAMMA_1, I_DELTA, I_NABLA_AD, I_OMEGA_ROT)
       is_defined = .TRUE.
    case default
       is_defined = .FALSE.
    end select

    ! Finish

    return

  end function is_defined

  !****

  function is_vacuum (this, pt)

    class(twopt_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    logical                          :: is_vacuum

    $ASSERT_DEBUG(pt%s >= this%s_i .AND. pt%s <= this%s_o,Invalid segment)

    ! Return whether the point is a vacuum

    is_vacuum = (1._WP - pt%x**2) == 0._WP

    ! Finish

    return

  end function is_vacuum

  !****

  function Delta_p (this, x_i, x_o)

    class(twopt_model_t), intent(in) :: this
    real(WP), intent(in)             :: x_i
    real(WP), intent(in)             :: x_o
    real(WP)                         :: Delta_p

    ! Evaluate the dimensionless p-mode frequency separation

    Delta_p = 0.5_WP/(SQRT(2._WP/this%Gamma_1)*(ASIN(x_o)-ASIN(x_i)))

    ! Finish

    return

  end function Delta_p

  !****

  function Delta_g (this, x_i, x_o, lambda)

    class(twopt_model_t), intent(in) :: this
    real(WP), intent(in)             :: x_i
    real(WP), intent(in)             :: x_o
    real(WP), intent(in)             :: lambda
    real(WP)                         :: Delta_g

    ! Evaluate the dimensionless g-mode inverse period separation

    Delta_g = 0._WP

    ! Finish

    return

  end function Delta_g

  !****

  function grid (this) result (gr)

    class(twopt_model_t), intent(in) :: this
    type(grid_t)                     :: gr

    ! Return the grid

    gr = this%gr

    ! Finish

    return

  end function grid

end module gyre_twopt_model
