! Module   : gyre_twopt_quintic_model
! Purpose  : stellar two-point model
!
! Copyright 2018-2019 The GYRE Team
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
     real(WP)     :: alpha(6)
     real(WP)     :: beta(4)
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
     procedure         :: dcoeff_U_
     procedure         :: dcoeff_c_1_
     procedure         :: P_ 
     procedure         :: rho_
     procedure         :: drho_dx_
     procedure         :: M_r_
     procedure         :: w_
     procedure         :: dw_dx_
     procedure         :: d2w_dx2_
     procedure         :: H3_
     procedure         :: dH3_dx_
     procedure         :: d2H3_dx2_
     procedure         :: H5_
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

    real(WP) :: M(4,6)

    ! Construct the twopt_quintic_model_t

    ml%gr = grid_t([0._WP, 1._WP])

    ml%Gamma_1 = ml_p%Gamma_1

    if (ml_p%uniform_rot) then
       ml%Omega_rot = uniform_Omega_rot(ml_p)
    else
       ml%Omega_rot = 0._WP
    endif

    ml%s_i = ml%gr%s_i()
    ml%s_o = ml%gr%s_o()

    ! Set up spline coefficients

    associate (P_c => ml_p%P_c, &
               P_s => ml_p%P_s, &
               rho_c => ml_p%rho_c, &
               rho_s => ml_p%rho_s)

      ! P-function

      ml%alpha(1) = P_c
      ml%alpha(2) = 0._WP
      ml%alpha(3) = -4._WP*PI*rho_c**2/3._WP
      ml%alpha(4) = 63._WP/(4._WP*PI) + 15._WP*(P_s - P_c - rho_c) + 2._WP*PI*rho_c**2/3._WP
      ml%alpha(5) = -rho_s
      ml%alpha(6) = P_s

      ! w-function

      M(1,1) = 0._WP
      M(1,2) = 0._WP
      M(1,3) = -84._WP
      M(1,4) = 0._WP
      M(1,5) = 0._WP
      M(1,6) = 0._WP

      M(2,1) = 2160._WP
      M(2,2) = 0._WP
      M(2,3) = 324._WP
      M(2,4) = -108._WP
      M(2,5) = 864._WP
      M(2,6) = -2160._WP

      M(3,1) = -360._WP
      M(3,2) = 0._WP
      M(3,3) = -12._WP
      M(3,4) = -24._WP
      M(3,5) = -144._WP
      M(3,6) = 360._WP

      M(4,1) = 60._WP
      M(4,2) = 0._WP
      M(4,3) = 2._WP
      M(4,4) = 4._WP
      M(4,5) = -60._WP
      M(4,6) = -60._WP

      ml%beta = PI/63._WP*MATMUL(M, ml%alpha)

    end associate

    ! Perform some sanity checks

    $ASSERT(ml%alpha(1) >= 0._WP,Negative central pressure)
    $ASSERT(ml%alpha(6) >= 0._WP,Negative surface pressure)

    $ASSERT(ml%beta(1) >= 0._WP,Negative central w)
    $ASSERT(ml%beta(4) >= 0._WP,Negative surface w)

    ! Finish

    return

  end function twopt_quintic_model_t_

  !****

  function coeff (this, i, pt)

    class(twopt_quintic_model_t), intent(in) :: this
    integer, intent(in)                      :: i
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: coeff

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
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: coeff

    real(WP) :: P
    real(WP) :: rho
    real(WP) :: w

    ! Evaluate the V_2 coefficient

    if (this%is_vacuum(pt)) then

       coeff = HUGE(0._WP)

    else

       associate (x => pt%x)

         P = this%P_(pt)
         rho = this%rho_(pt)
         w = this%w_(pt)
         
         coeff = SQRT(w)*rho/P

       end associate

    end if

    ! Finish

    return

  end function coeff_V_2_

  !****

  function coeff_As_ (this, pt) result (coeff)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: coeff

    real(WP) :: P
    real(WP) :: rho
    real(WP) :: drho_dx
    real(WP) :: w

    ! Evaluate the As coefficient

    if (this%is_vacuum(pt)) then

       coeff = HUGE(0._WP)

    else

       associate (x => pt%x, &
                  Gamma_1 => this%Gamma_1)

         P = this%P_(pt)
         rho = this%rho_(pt)
         drho_dx = this%drho_dx_(pt)
         w = this%w_(pt)
         
         coeff = -x**2*SQRT(w)*rho/(Gamma_1*P) - x*drho_dx/rho

       end associate

    end if

    ! Finish

    return

  end function coeff_As_

  !****

  function coeff_U_ (this, pt) result (coeff)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: coeff

    real(WP) :: rho
    real(WP) :: w

    ! Evaluate the U coefficient

    associate (x => pt%x)

      rho = this%rho_(pt)
      w = this%w_(pt)

      coeff = 4._WP*PI*rho/SQRT(w)

    end associate

    ! Finish

    return

  end function coeff_U_

  !****

  function coeff_c_1_ (this, pt) result (coeff)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: coeff

    real(WP) :: w

    ! Evaluate the c_1 coefficient

    associate (x => pt%x)

      w = this%w_(pt)

      coeff = 1._WP/w

    end associate

    ! Finish

    return

  end function coeff_c_1_

  !****

  function dcoeff (this, i, pt)

    class(twopt_quintic_model_t), intent(in) :: this
    integer, intent(in)                      :: i
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: dcoeff

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

    real(WP) :: rho
    real(WP) :: drho_dx

    ! Evaluate the logarithmic derivative of the V_2 coefficient

    associate (x => pt%x)

      rho = this%rho_(pt)
      drho_dx = this%drho_dx_(pt)
      
      dcoeff = this%coeff_U_(pt) + x*drho_dx/rho + this%coeff_V_2_(pt)*x**2 - 3._WP

    end associate

    ! Finish

    return

  end function dcoeff_V_2_

  !****

  function dcoeff_U_ (this, pt) result (dcoeff)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: dcoeff

    real(WP) :: rho
    real(WP) :: drho_dx
    real(WP) :: w
    real(WP) :: dw_dx

    ! Evaluate the logarithmic derivative of the U coefficient

    associate (x => pt%x)

      rho = this%rho_(pt)
      drho_dx = this%drho_dx_(pt)

      w = this%w_(pt)
      dw_dx = this%dw_dx_(pt)

      dcoeff = x*drho_dx/rho - x*dw_dx/(2._WP*w)

    end associate


    ! Finish

    return

  end function dcoeff_U_

  !****

  function dcoeff_c_1_ (this, pt) result (dcoeff)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: dcoeff

    ! Evaluate the logarithmic derivative of the c_1 coefficient

    dcoeff = 3._WP - this%coeff_U_(pt)

    ! Finish

    return

  end function dcoeff_c_1_

  !****

  function P_ (this, pt) result (P)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: P

    real(WP) :: H5(6)

    ! Evaluate the dimensionless pressure

    H5 = this%H5_(pt)

    P = DOT_PRODUCT(this%alpha, H5)

    ! Finish

    return

  end function P_

  !****

  function rho_ (this, pt) result (rho)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: rho

    real(WP) :: w
    real(WP) :: dw_dx

    ! Evaluate the dimensionless density

    associate (x => pt%x)

      w = this%w_(pt)
      dw_dx = this%dw_dx_(pt)

      rho = 3._WP*SQRT(w)/(4._WP*PI)*(1._WP + x*dw_dx/(6._WP*w))

    end associate

    ! Finish

    return

  end function rho_

  !****

  function drho_dx_ (this, pt) result (drho_dx)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: drho_dx

    real(WP) :: w
    real(WP) :: dw_dx
    real(WP) :: d2w_dx2

    ! Evaluate the derivative of the dimensionless density

    associate (x => pt%x)

      w = this%w_(pt)
      dw_dx = this%dw_dx_(pt)
      d2w_dx2 = this%d2w_dx2_(pt)

      drho_dx = (2._WP*w*(4._WP*dw_dx + x*d2w_dx2) - x*dw_dx**2)/(16._WP*PI*w*SQRT(w))

    end associate

    ! Finish

    return

  end function drho_dx_

  !****

  function M_r_ (this, pt) result (M_r)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: M_r

    real(WP) :: w

    ! Evaluate the dimensionless mass coordinate

    associate (x => pt%x)

      w = this%w_(pt)
         
      M_r = x**3*SQRT(w)

    end associate

    ! Finish

    return

  end function M_r_

  !****

  function w_ (this, pt) result (w)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: w

    real(WP) :: H3(4)

    ! Evaluate the w = x**6 M_r**2 function

    H3 = this%H3_(pt)

    w = DOT_PRODUCT(this%beta, H3)

    if (w < 0._WP) then
       write(*,*) 'x, w:', pt%x, w
       $ABORT(Negative w)
    endif

    ! Finish

    return

  end function w_

  !****

  function dw_dx_ (this, pt) result (dw_dx)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: dw_dx

    real(WP) :: dH3_dx(4)

    ! Evaluate the derivative of the w function

    dH3_dx = this%dH3_dx_(pt)

    dw_dx = DOT_PRODUCT(this%beta, dH3_dx)

    ! Finish

    return

  end function dw_dx_

  !****

  function d2w_dx2_ (this, pt) result (d2w_dx2)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: d2w_dx2

    real(WP) :: d2H3_dx2(4)

    ! Evaluate the second derivative of the w function

    d2H3_dx2 = this%d2H3_dx2_(pt)

    d2w_dx2 = DOT_PRODUCT(this%beta, d2H3_dx2)

    ! Finish

    return
 
  end function d2w_dx2_

  !****

  function H3_ (this, pt) result (H3)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: H3(4)

    ! Evaluate the cubic Hermite basis functions

    associate (x => pt%x)
    
      H3(1) = 1._WP - 3._WP*x**2 + 2._WP*x**3
      H3(2) = x - 2._WP*x**2 + x**3
      H3(3) = -x**2 + x**3
      H3(4) = 3._WP*x**2 - 2._WP*x**3

    end associate

    ! Finish

    return

  end function H3_

  !****

  function dH3_dx_ (this, pt) result (dH3_dx)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: dH3_dx(4)

    ! Evaluate the derivative of the cubic Hermite basis functions

    associate (x => pt%x)
    
      dH3_dx(1) = -6._WP*x + 6._WP*x**2
      dH3_dx(2) = 1._WP - 4._WP*x + 3._WP*x**2
      dH3_dx(3) = -2._WP*x + 3._WP*x**2
      dH3_dx(4) = 6._WP*x - 6._WP*x**2

    end associate

    ! Finish

    return

  end function dH3_dx_

  !****

  function d2H3_dx2_ (this, pt) result (d2H3_dx2)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: d2H3_dx2(4)

    ! Evaluate the second derivative of the cubic Hermite basis functions

    associate (x => pt%x)
    
      d2H3_dx2(1) = -6._WP + 12._WP*x
      d2H3_dx2(2) = -4._WP + 6._WP*x
      d2H3_dx2(3) = -2._WP + 6._WP*x
      d2H3_dx2(4) = 6._WP - 12._WP*x

    end associate

    ! Finish

    return

  end function d2H3_dx2_

  !****

  function H5_ (this, pt) result (H5)

    class(twopt_quintic_model_t), intent(in) :: this
    type(point_t), intent(in)                :: pt
    real(WP)                                 :: H5(6)

    ! Evaluate the quintic Hermite basis functions

    associate (x => pt%x)
    
      H5(1) = 1._WP - 10._WP*x**3 + 15._WP*x**4 - 6._WP*x**5
      H5(2) = x - 6._WP*x**3 + 8._WP*x**4 - 3._WP*x**5
      H5(3) = x**2/2._WP - 3._WP*x**3/2._WP + 3._WP*x**4/2._WP - x**5/2._WP
      H5(4) = x**3/2._WP - x**4 + x**5/2._WP
      H5(5) = -4._WP*x**3 + 7._WP*x**4 - 3._WP*x**5
      H5(6) = 10._WP*x**3 - 15._WP*x**4 + 6._WP*x**5
    
    end associate

    ! Finish

    return

  end function H5_

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

    is_vacuum = this%P_(pt) == 0._WP

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
