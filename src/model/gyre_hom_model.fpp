! Module   : gyre_hom_model
! Purpose  : stellar homogeneous compressible model
!
! Copyright 2013-2017 Rich Townsend
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

module gyre_hom_model

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

  type, extends (model_t) :: hom_model_t
     private
     type(grid_t) :: gr
     real(WP)     :: Gamma_1
     real(WP)     :: Omega_rot
     integer      :: s_i
     integer      :: s_o
   contains
     private
     procedure, public :: coeff
     procedure         :: coeff_V_2_
     procedure         :: coeff_As_
     procedure, public :: dcoeff
     procedure         :: dcoeff_V_2_
     procedure         :: dcoeff_As_
     procedure, public :: is_defined
     procedure, public :: is_vacuum
     procedure, public :: Delta_p
     procedure, public :: Delta_g
     procedure, public :: grid
  end type hom_model_t

  ! Interfaces

  interface hom_model_t
     module procedure hom_model_t_
  end interface hom_model_t

  ! Access specifiers

  private

  public :: hom_model_t

  ! Procedures

contains

  function hom_model_t_ (ml_p) result (ml)

    use gyre_grid_weights
  
    type(model_par_t), intent(in) :: ml_p
    type(hom_model_t)             :: ml

    real(WP), allocatable :: w(:)

    ! Construct the hom_model_t

    select case (ml_p%grid_type)
    case ('UNI')
       w = uni_weights(ml_p%n)
    case ('GEO')
       w = geo_weights(ml_p%n, ml_p%s)
    case ('LOG')
       w = log_weights(ml_p%n, ml_p%s)
    case default
       $ABORT(Invalid grid_type)
    end select

    ml%gr = grid_t(w, ml_p%x_i, ml_p%x_o)

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

  end function hom_model_t_

  !****

  function coeff (this, i, pt)

    class(hom_model_t), intent(in) :: this
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
       coeff = 3._WP
    case (I_C_1)
       coeff = 1._WP
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

    class(hom_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: coeff

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),V_2 evaluation at vacuum point)

    ! Evaluate the V_2 coefficient

    coeff = 2._WP/(1._WP - pt%x**2)

    ! Finish

    return

  end function coeff_V_2_

  !****

  function coeff_As_ (this, pt) result (coeff)

    class(hom_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: coeff

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),As evaluation at vacuum point)

    ! Evaluate the As coefficient

    coeff = -this%coeff_V_2_(pt)*pt%x**2/this%Gamma_1

    ! Finish

    return

  end function coeff_As_

  !****

  function dcoeff (this, i, pt)

    class(hom_model_t), intent(in) :: this
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
       dcoeff = 0._WP
    case (I_C_1)
       dcoeff = 0._WP
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

    class(hom_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: dcoeff

    ! Evaluate the logarithmic derivative of the V_2 coefficient

    dcoeff = 2.*WP*pt%x**2/(1._WP - pt%x**2)

    ! Finish

    return

  end function dcoeff_V_2_

  !****

  function dcoeff_As_ (this, pt) result (dcoeff)

    class(hom_model_t), intent(in) :: this
    type(point_t), intent(in)      :: pt
    real(WP)                       :: dcoeff

    ! Evaluate the logarithmic derivative of the As coefficient
    
    dcoeff = -this%dcoeff_V_2_(pt)*pt%x**2/this%Gamma_1

    ! Finish

    return

  end function dcoeff_As_

  !****

  function is_defined (this, i)

    class(hom_model_t), intent(in) :: this
    integer, intent(in)            :: i
    logical                        :: is_defined

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

    class(hom_model_t), intent(in) :: this
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

    class(hom_model_t), intent(in) :: this
    real(WP), intent(in)           :: x_i
    real(WP), intent(in)           :: x_o
    real(WP)                       :: Delta_p

    ! Evaluate the dimensionless p-mode frequency separation

    Delta_p = 0.5_WP/(SQRT(2._WP/this%Gamma_1)*(ASIN(x_o)-ASIN(x_i)))

    ! Finish

    return

  end function Delta_p

  !****

  function Delta_g (this, x_i, x_o, lambda)

    class(hom_model_t), intent(in) :: this
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

    class(hom_model_t), intent(in) :: this
    type(grid_t)                   :: gr

    ! Return the grid

    gr = this%gr

    ! Finish

    return

  end function grid

end module gyre_hom_model
