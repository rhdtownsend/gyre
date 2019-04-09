! Module   : gyre_twopt_quintic_model
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

module gyre_twopt_quintic_model

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

  type, extends (model_t) :: twopt_quintic_model_t
     private
     type(grid_t) :: gr
     real(WP)     :: P_cent_2pt
     real(WP)     :: P_surf_2pt
     real(WP)     :: rho_cent_2pt
     real(WP)     :: rho_surf_2pt
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
     procedure         :: drho_
     procedure         :: m_r_
     procedure         :: P_
     procedure         :: dcoeff_U_
     procedure         :: dcoeff_c_1_
     procedure         :: alpha_
     procedure         :: beta_
     procedure         :: hermite_five
     procedure         :: hermite_three
     procedure         :: w_
     procedure         :: dw_
     procedure         :: ddw_
     procedure, public :: is_defined
     procedure, public :: is_vacuum
     procedure, public :: Delta_p
     procedure, public :: Delta_g
     procedure, public :: grid
  end type twopt_quintic_model_t

  ! Interfaces

  interface twopt_quintic_model_t
     module procedure twopt_quintic_model_t_
  end interface twopt_quintic_model_t

  
  ! Access specifiers

  private

  public :: twopt_quintic_model_t

  ! Procedures

contains

  function twopt_quintic_model_t_ (ml_p) result (ml)

    use gyre_grid_weights
  
    type(model_par_t), intent(in) :: ml_p
    type(twopt_quintic_model_t)   :: ml


    ! Construct the twopt_quintic_model_t

    ml%gr = grid_t([0._WP, 1._WP])

    ml%Gamma_1 = ml_p%Gamma_1

    ml%P_cent_2pt = ml_p%P_cent_2pt
    ml%P_surf_2pt = ml_p%P_surf_2pt
    ml%rho_cent_2pt = ml_p%rho_cent_2pt
    ml%rho_surf_2pt = ml_p%rho_surf_2pt

    if (ml_p%uniform_rot) then
       ml%Omega_rot = uniform_Omega_rot(ml_p)
    else
       ml%Omega_rot = 0._WP
    endif

    ml%s_i = ml%gr%s_i()
    ml%s_o = ml%gr%s_o()

    ! Finish

    return

  end function twopt_quintic_model_t_

  !****

  function coeff (this, i, pt)

    class(twopt_quintic_model_t), intent(in) :: this
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

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: coeff

    ! Evaluate the V_2 coefficient

    associate ( &
         x => pt%x, &
         w => this%w_(pt))

    IF (is_vacuum(this,pt) .eqv. .FALSE.) THEN
        coeff = (SQRT(w) * this%rho_(pt)) / (this%P_(pt))
    ELSE IF (is_vacuum(this,pt) .eqv. .TRUE.) THEN
        coeff = HUGE(0._WP)
    END IF

    end associate

    ! Finish

    return

  end function coeff_V_2_

  !****

  function coeff_As_ (this, pt) result (coeff)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: coeff

    ! Evaluate the As coefficient

    coeff = 3._WP - this%coeff_U_(pt) - this%coeff_V_2_(pt)*pt%x**2/this%Gamma_1 - this%dcoeff_U_(pt)

    ! Finish

    return

  end function coeff_As_

  !****

  function coeff_U_ (this, pt) result (coeff)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: coeff

    ! Evaluate the U coefficient

    associate ( &
         x => pt%x, &
         w => this%w_(pt))

      coeff = (4*PI *  this%rho_(pt))/SQRT(w)

    end associate

    ! Finish

    return

  end function coeff_U_

  !****

  function coeff_c_1_ (this, pt) result (coeff)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: coeff

    ! Evaluate the c_1 coefficient

    associate ( &
         x => pt%x, &
         w => this%w_(pt))

      coeff = 1._WP/SQRT(w)

    end associate

    ! Finish

    return

  end function coeff_c_1_

  !****

  function dcoeff (this, i, pt)

    class(twopt_quintic_model_t), intent(in) :: this
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

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: dcoeff

    ! Evaluate the logarithmic derivative of the V_2 coefficient

    associate ( &
         x => pt%x)

      
      dcoeff = this%coeff_U_(pt) + (x/this%rho_(pt))*this%drho_(pt) + &
        this%coeff_V_2_(pt)*x**2 - 3._WP

    end associate

    ! Finish

    return

  end function dcoeff_V_2_

  !****

  function rho_ (this, pt) result (rho)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: rho

    ! Evaluate the polynomial for density function

    associate ( &
         x => pt%x, &
         w => this%w_(pt), &
         dw => this%dw_(pt))

      rho = ((3._WP*SQRT(w))/(4*PI)) * (1._WP + (x*dw)/(6._WP*w))

    end associate

  end function rho_

  !****

  function drho_ (this, pt) result (drho)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: drho

    ! Evaluate the polynomial for density function

    associate ( &
         x => pt%x, &
         w => this%w_(pt), &
         dw => this%dw_(pt), &
         ddw => this%ddw_(pt))

      drho = (6._WP*w + dw*(2*w + x - 2._WP*x*dw) + 2._WP*x*w*ddw)/(16._WP*w**(3._WP/2._WP)*PI)

    end associate

  end function drho_

  !****

  function m_r_ (this, pt) result (m_r)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: m_r

    ! Evaluate the polynomial for m_n profile

    associate ( &
         x => pt%x, &
         w => this%w_(pt))
         
      m_r = x**3 * SQRT(w)

    end associate

  end function m_r_

  !****

  function P_ (this, pt) result (Pressure)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: Pressure

    ! Evaluate the polynomial for Pressure profile  

    associate ( &
         x => pt%x)

      Pressure = (this%alpha_(1, pt)*this%hermite_five(1,pt)&
            + this%alpha_(2, pt)*this%hermite_five(2,pt)&
            + this%alpha_(3, pt)*this%hermite_five(3,pt)&
            + this%alpha_(4, pt)*this%hermite_five(4,pt)&
            + this%alpha_(5, pt)*this%hermite_five(5,pt)&
            + this%alpha_(6, pt)*this%hermite_five(6,pt))

    end associate

  end function P_

  !****

  function dcoeff_U_ (this, pt) result (dcoeff)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: dcoeff

    ! Evaluate the logarithmic derivative of the U coefficient

    associate ( &
         x => pt%x, &
         rho => this%rho_(pt), &
         drho => this%drho_(pt), &
         w => this%w_(pt), &
         dw => this%dw_(pt))

      dcoeff = -(x*dw)/(2._WP*w) + (x*drho)/rho

    end associate


    ! Finish

    return

  end function dcoeff_U_

  !****

  function dcoeff_c_1_ (this, pt) result (dcoeff)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: dcoeff

    ! Evaluate the logarithmic derivative of the c_1 coefficient

    dcoeff = 3._WP - this%coeff_U_(pt)

    ! Finish

    return

  end function dcoeff_c_1_

  !****

  function alpha_ (this, n, pt) result (alpha)

    class(twopt_quintic_model_t), intent(in) :: this
    integer, intent(in)              :: n
    type(point_t), intent(in)        :: pt
    real(WP)                         :: alpha

    ! Evaluate alpha coefficients

    associate (x => pt%x)

    select case(n)
    case(1)
        alpha = this%P_cent_2pt
    case(2)
        alpha = 0._WP
    case(3)
        alpha = -(4._WP * PI / 3._WP) * (this%rho_cent_2pt)**2
    case(4)
        alpha = (63._WP/(4._WP * PI)) + (2._WP * pi/3._WP)*(this%rho_cent_2pt)**2 + &
            15._WP*(this%P_surf_2pt - this%P_cent_2pt - this%rho_surf_2pt)
    case(5)
        alpha = -this%rho_surf_2pt
    case(6)
        alpha = this%P_surf_2pt
    end select

    end associate
    
    ! Finish

    return

  end function alpha_

  !****

  function beta_ (this, n, pt) result (beta)

    class(twopt_quintic_model_t), intent(in) :: this
    integer, intent(in)              :: n
    type(point_t), intent(in)        :: pt
    real(WP)                         :: beta
    real(WP)                         :: beta_sub

    ! Evaluate beta coefficients

    associate (x => pt%x)

    select case(n)
    case(1)
        beta_sub = -84._WP * this%alpha_(3, pt)
    case(2)
        beta_sub = 2160._WP*this%alpha_(1, pt) + 324._WP*this%alpha_(3, pt) - &
            108._WP*this%alpha_(4, pt) + 864._WP*this%alpha_(5, pt) - 2160._WP*this%alpha_(6, pt)
    case(3)
        beta_sub = -360._WP*this%alpha_(1, pt) - 12._WP*this%alpha_(3, pt) - &
            12._WP*this%alpha_(4, pt) - 144._WP*this%alpha_(5, pt) + 360._WP*this%alpha_(6, pt)
    case(4)
        beta_sub = 60._WP*this%alpha_(1, pt) + 2._WP*this%alpha_(3, pt) + &
            4._WP*this%alpha_(4, pt) - 60._WP*this%alpha_(5, pt) - 60._WP*this%alpha_(6, pt)
    end select

    beta = (PI/63._WP) * beta_sub

    end associate
    
    ! Finish

    return

  end function beta_

  !****

  function w_(this, pt) result (w)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: w

    associate (x => pt%x)

    w = this%beta_(1,pt)*this%hermite_three(1,pt) + &
      this%beta_(2,pt)*this%hermite_three(2,pt) + &
      this%beta_(3,pt)*this%hermite_three(3,pt) + &
      this%beta_(4,pt)*this%hermite_three(4,pt)

    end associate

    return

  end function w_

  !****

  function dw_(this, pt) result (dw)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: dw

    associate ( &
         x => pt%x, &
         a1 => this%alpha_(1,pt), &
         a2 => this%alpha_(2,pt), &
         a3 => this%alpha_(3,pt), &
         a4 => this%alpha_(4,pt), &
         a5 => this%alpha_(5,pt), &
         a6 => this%alpha_(6,pt))

    dw = (4._WP * PI / 21._WP) * (9._WP*(20._WP*a1 + 3._WP*a3 - a4 + 8._WP*a5 - 20._WP*a6) + &
        35._WP* x**2 *(12._WP*a1 + a3 - a4 + 6._WP*a5 - 12._WP*a6) - 21._WP*x*(30._WP*a1 + &
        3._WP*a3 - 2._WP*a4 + 14._WP*a5 - 30._WP*a6))

    end associate

    return

  end function dw_

  !****

  function ddw_(this, pt) result (ddw)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)        :: pt
    real(WP)                         :: ddw

    associate ( &
         x => pt%x, &
         a1 => this%alpha_(1,pt), &
         a2 => this%alpha_(2,pt), &
         a3 => this%alpha_(3,pt), &
         a4 => this%alpha_(4,pt), &
         a5 => this%alpha_(5,pt), &
         a6 => this%alpha_(6,pt))

    ddw = (40._WP * PI * x/3._WP)*(12._WP*a1 + a3 - a4 + 6._WP*a5 - 12._WP*a6) - &
        4._WP*PI*(30._WP*a1 + 3._WP*a3 - 2._WP*(a4 - 7._WP*a5 + 15._WP*a6))

    end associate

    return

  end function ddw_


  !****

  function hermite_five (this, n, pt) result (basis)

    class(twopt_quintic_model_t), intent(in) :: this
    integer, intent(in)              :: n
    type(point_t), intent(in)        :: pt
    real(WP)                         :: basis


    ! Evaluate the n'th Quintic Hermite Basis Function

    associate (x => pt%x)
    
    select case (n)
    case (1)
       basis = 1._WP - 10._WP*x**3 + 15._WP*x**4 - 6._WP*x**5
    case (2)
       basis = x - 6._WP*x**3 + 8._WP*x**4 - 3._WP*x**5
    case (3)
       basis = (x**2)/2._WP - 3._WP*(x**3)/2._WP + 3._WP*(x**4)/2._WP - (x**5)/2._WP
    case (4)
       basis = (x**3)/2._WP - x**4 + (x**5)/2
    case (5)
       basis = -4._WP*x**3 + 7._WP*(x**4) - 3._WP*x**5
    case (6)
       basis = 10._WP*x**3 - 15._WP*x**4 + 6._WP*x**5
    
    end select

    end associate
    ! Finish

    return

  end function hermite_five

  !****

  function hermite_three (this, n, pt) result (basis)

    class(twopt_quintic_model_t), intent(in) :: this
    integer, intent(in)              :: n
    type(point_t), intent(in)        :: pt
    real(WP)                         :: basis


    ! Evaluate the n'th Cubic Hermite Basis Function

    associate (x => pt%x)
    
    select case (n)
    case (1)
       basis = 1._WP - 3._WP*x**2 + 2._WP*x**3
    case (2)
       basis = x - 2._WP*x**2 + x**3
    case (3)
       basis = -x**2 + x**3
    case (4)
       basis = 3._WP*x**2 - 2._WP*x**3
    
    end select

    end associate
    ! Finish

    return

  end function hermite_three

  !****

  function is_defined (this, i)

    class(twopt_quintic_model_t), intent(in) :: this
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

    class(twopt_quintic_model_t), intent(in) :: this
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

    class(twopt_quintic_model_t), intent(in) :: this
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

    class(twopt_quintic_model_t), intent(in) :: this
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

    class(twopt_quintic_model_t), intent(in) :: this
    type(grid_t)                     :: gr

    ! Return the grid

    gr = this%gr

    ! Finish

    return

  end function grid

end module gyre_twopt_quintic_model
