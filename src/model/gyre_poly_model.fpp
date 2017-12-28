! Module   : gyre_poly_model
! Purpose  : stellar polytropic model
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

module gyre_poly_model

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_grid
  use gyre_interp
  use gyre_model
  use gyre_point

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (model_t) :: poly_model_t
     private
     type(grid_t)                  :: gr
     type(r_interp_t), allocatable :: in_Theta(:)
     type(r_interp_t), allocatable :: in_dTheta(:)
     real(WP), allocatable         :: n_poly(:)
     real(WP), allocatable         :: mu_i(:)
     real(WP), allocatable         :: v_i(:)
     real(WP), allocatable         :: t(:)
     real(WP), allocatable         :: B(:)
     real(WP)                      :: mu_s
     real(WP)                      :: xi_s
     real(WP)                      :: Gamma_1
     real(WP)                      :: Omega_rot
     integer                       :: s_i
     integer                       :: s_o
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
     procedure         :: mu_
     procedure, public :: is_defined
     procedure, public :: is_vacuum
     procedure, public :: Delta_p
     procedure, public :: Delta_g
     procedure, public :: grid
  end type poly_model_t

  ! Interfaces

  interface poly_model_t
     module procedure poly_model_t_
  end interface poly_model_t

  ! Access specifiers

  private

  public :: poly_model_t

  ! Procedures

contains

  function poly_model_t_ (xi, Theta, dTheta, n_poly, Delta_d, Gamma_1, Omega_rot) result (ml)

    real(WP), intent(in) :: xi(:)
    real(WP), intent(in) :: Theta(:)
    real(WP), intent(in) :: dTheta(:)
    real(WP), intent(in) :: n_poly(:)
    real(WP), intent(in) :: Delta_d(:)
    real(WP), intent(in) :: Gamma_1
    real(WP), intent(in) :: Omega_rot
    type(poly_model_t)   :: ml

    integer  :: n_k
    real(WP) :: x(SIZE(xi))
    integer  :: k_i
    integer  :: k_o
    integer  :: k_i_prev
    integer  :: k_o_prev
    integer  :: s
    integer  :: i
    real(WP) :: v_o_prev
    real(WP) :: d2Theta(SIZE(xi))

    $CHECK_BOUNDS(SIZE(Theta),SIZE(xi))
    $CHECK_BOUNDS(SIZE(dTheta),SIZE(xi))

    $CHECK_BOUNDS(SIZE(Delta_d),SIZE(n_poly)-1)

    ! Construct the poly_model_t from the Lane-Emden solutions Theta,
    ! dTheta/dxi. Per-segment polytropic indices and density jumps are
    ! supplied in n_poly and Delta_d, respectively

    ! Create the grid

    n_k = SIZE(xi)

    ml%xi_s = xi(n_k)

    x = xi/ml%xi_s

    ml%gr = grid_t(x)

    ml%s_i = ml%gr%s_i()
    ml%s_o = ml%gr%s_o()

    $CHECK_BOUNDS(SIZE(n_poly),ml%s_o-ml%s_i+1)

    ! Allocate arrays

    allocate(ml%in_Theta(ml%s_i:ml%s_o))
    allocate(ml%in_dTheta(ml%s_i:ml%s_o))

    allocate(ml%n_poly(ml%s_i:ml%s_o))
    allocate(ml%mu_i(ml%s_i:ml%s_o))
    allocate(ml%v_i(ml%s_i:ml%s_o))
    allocate(ml%t(ml%s_i:ml%s_o))
    allocate(ml%B(ml%s_i:ml%s_o))

    ml%n_poly = n_poly

    ! Set up per-segment mu_i, v_i, t and B data

    ml%mu_i(ml%s_i) = 0._WP
    ml%v_i(ml%s_i) = 0._WP
    ml%t(ml%s_i) = 1._WP
    ml%B(ml%s_i) = 1._WP

    k_i = ml%gr%k_s_i(ml%s_i)
    k_o = ml%gr%k_s_o(ml%s_i)

    seg_data_loop : do s = ml%s_i+1, ml%s_o

       k_i_prev = k_i
       k_o_prev = k_o

       k_i = ml%gr%k_s_i(s)
       k_o = ml%gr%k_s_o(s)

       i = s - ml%s_i + 1

       v_o_prev = xi(k_o_prev)**2*dTheta(k_o_prev)

       ml%mu_i(s) = ml%mu_i(s-1) - (v_o_prev - ml%v_i(s-1))*ml%t(s-1)/ml%B(s-1)

       ml%t(s) = ml%t(s-1)*EXP(ml%n_poly(s-1)*LOG(Theta(k_o_prev)) + Delta_d(i-1))

       ml%v_i(s) = xi(k_i)**2*dTheta(k_i)

       ml%B(s) = (dTheta(k_i)/dTheta(k_o_prev))*(ml%t(s)/ml%t(s-1))*ml%B(s-1)

    end do seg_data_loop

    v_o_prev = xi(k_o)**2*dTheta(k_o)

    ml%mu_s = ml%mu_i(s-1) - (v_o_prev - ml%v_i(s-1))*ml%t(s-1)/ml%B(s-1)

    ! Set up per-segment splines

    seg_spline_loop : do s = ml%s_i, ml%s_o

       k_i = ml%gr%k_s_i(s)
       k_o = ml%gr%k_s_o(s)

       if (ml%n_poly(s) /= 0._WP) then

          where (xi(k_i:k_o) /= 0._WP)
             d2Theta(k_i:k_o) = -2._WP*dTheta(k_i:k_o)/xi(k_i:k_o) - ml%B(s)*Theta(k_i:k_o)**ml%n_poly(s)
          elsewhere
             d2Theta(k_i:k_o) = -1._WP/3._WP
          end where

       else

          d2Theta(k_i:k_o) = -1._WP/3._WP

       endif

       ml%in_Theta(s) = r_interp_t(x(k_i:k_o), Theta(k_i:k_o), dTheta(k_i:k_o)*ml%xi_s)
       ml%in_dTheta(s) = r_interp_t(x(k_i:k_o), dTheta(k_i:k_o), d2Theta(k_i:k_o)*ml%xi_s)

    end do seg_spline_loop

    ! Other initializations

    ml%Gamma_1 = Gamma_1
    ml%Omega_rot = Omega_rot

    ! Finish

    return

  end function poly_model_t_

  !****

  function coeff (this, i, pt)

    class(poly_model_t), intent(in) :: this
    integer, intent(in)             :: i
    type(point_t), intent(in)       :: pt
    real(WP)                        :: coeff

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

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: coeff

    real(WP) :: xi
    real(WP) :: Theta
    real(WP) :: dTheta

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),V_2 evaluation at vacuum point)

    ! Evaluate the V_2 coefficient

    if (pt%x /= 0._WP) then

       xi = pt%x*this%xi_s

       Theta = this%in_Theta(pt%s)%f(pt%x)
       dTheta = this%in_dTheta(pt%s)%f(pt%x)

       coeff = -(this%n_poly(pt%s) + 1._WP)*this%xi_s**2*dTheta/(Theta*xi)

    else

       coeff = (this%n_poly(pt%s) + 1._WP)*this%xi_s**2/3._WP

    endif

    ! Finish

    return

  end function coeff_V_2_

  !****

  function coeff_As_ (this, pt) result (coeff)

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: coeff

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),As evaluation at vacuum point)

    ! Evaluate the As coefficient

    coeff = this%coeff_V_2_(pt)*pt%x**2 * &
         (this%n_poly(pt%s)/(this%n_poly(pt%s) + 1._WP) - 1._WP/this%Gamma_1)

    ! Finish

    return

  end function coeff_As_

  !****

  function coeff_U_ (this, pt) result (coeff)

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: coeff

    real(WP) :: xi
    real(WP) :: Theta

    ! Evaluate the U coefficient

    if (pt%x /= 0._WP) then

       xi = pt%x*this%xi_s

       Theta = this%in_Theta(pt%s)%f(pt%x)

       coeff = xi**3*this%t(pt%s)*Theta**this%n_poly(pt%s)/this%mu_(pt)

    else

       coeff = 3._WP

    endif

    ! Finish

    return

  end function coeff_U_

  !****

  function coeff_c_1_ (this, pt) result (coeff)

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: coeff

    ! Evaluate the c_1 coefficient

    if (pt%x /= 0._WP) then

       coeff = pt%x**3/(this%mu_(pt)/this%mu_s)

    else

       coeff = 3._WP*this%mu_s/this%xi_s**3

    endif

    ! Finish

    return

  end function coeff_c_1_

  !****

  function dcoeff (this, i, pt)

    class(poly_model_t), intent(in) :: this
    integer, intent(in)             :: i
    type(point_t), intent(in)       :: pt
    real(WP)                        :: dcoeff

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

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: dcoeff

    real(WP) :: xi
    real(WP) :: Theta
    real(WP) :: dTheta

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),dV_2 evaluation at vacuum point)

    ! Evaluate the logarithmic derivative of the V_2 coefficient

    if (pt%x /= 0._WP) then

       xi = pt%x*this%xi_s

       Theta = this%in_Theta(pt%s)%f(pt%x)
       dTheta = this%in_dTheta(pt%s)%f(pt%x)

       dcoeff = -3._WP - xi*dTheta/Theta - this%B(pt%s)*xi*Theta**this%n_poly(pt%s)/dTheta

    else

       dcoeff = 0._WP

    endif

    ! Finish

    return

  end function dcoeff_V_2_

  !****

  function dcoeff_As_ (this, pt) result (dcoeff)

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: dcoeff

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),dAs evaluation at vacuum point)

    ! Evaluate the logarithmic derivative of the As coefficient
    
    dcoeff = this%dcoeff_V_2_(pt)*pt%x**2 * &
         (this%n_poly(pt%s)/(this%n_poly(pt%s) + 1._WP) - 1._WP/this%Gamma_1)

    ! Finish

    return

  end function dcoeff_As_

  !****

  function dcoeff_U_ (this, pt) result (dcoeff)

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: dcoeff

    real(WP) :: xi
    real(WP) :: Theta
    real(WP) :: dTheta

    $ASSERT_DEBUG(.NOT. this%is_vacuum(pt),dU evaluation at vacuum point)

    ! Evaluate the logarithmic derivative of the U coefficient
    
    if (pt%x /= 0._WP) then

       xi = pt%x*this%xi_s

       Theta = this%in_Theta(pt%s)%f(pt%x)
       dTheta = this%in_dTheta(pt%s)%f(pt%x)
       
       dcoeff = 3._WP + this%n_poly(pt%s)*xi*dTheta/Theta - this%coeff_U_(pt)

    else

       dcoeff = 0._WP

    endif

    ! Finish

    return

  end function dcoeff_U_

  !****

  function dcoeff_c_1_ (this, pt) result (dcoeff)

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: dcoeff

    ! Evaluate the logarithmic derivative of the c_1 coefficient
    
    dcoeff = 3._WP - this%coeff_U_(pt)

    ! Finish

    return

  end function dcoeff_c_1_

  !****

  function mu_ (this, pt) result (mu)

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: mu

    real(WP) :: v

    ! Evaluate the mass coordinate mu

    v = (this%xi_s*pt%x)**2*this%in_dTheta(pt%s)%f(pt%x)

    mu = this%mu_i(pt%s) - (v - this%v_i(pt%s))*this%t(pt%s)/this%B(pt%s)

    ! Finish

    return

  end function mu_

  !****

  function is_defined (this, i)

    class(poly_model_t), intent(in) :: this
    integer, intent(in)             :: i
    logical                         :: is_defined

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

    class(poly_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    logical                         :: is_vacuum

    $ASSERT_DEBUG(pt%s >= this%s_i .AND. pt%s <= this%s_o,Invalid segment)

    ! Return whether the point is a vacuum

    is_vacuum = this%in_Theta(pt%s)%f(pt%x) == 0._WP

    ! Finish

    return

  end function is_vacuum

  !****

  function Delta_p (this, x_i, x_o)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x_i
    real(WP), intent(in)            :: x_o
    real(WP)                        :: Delta_p

    type(grid_t)  :: gr
    real(WP)      :: I
    integer       :: s
    type(point_t) :: pt
    integer       :: k_i
    integer       :: k_o
    integer       :: k
    real(WP)      :: V_2
    real(WP)      :: c_1
    real(WP)      :: Gamma_1

    ! Evaluate the dimensionless g-mode inverse period separation,
    ! using a midpoint quadrature rule since the integrand typically
    ! diverges at the surface

    ! First, create the nested grid

    gr = grid_t(this%gr, x_i, x_o)

    ! Now evaluate the integrand segment by segment

    I = 0._WP

    seg_loop : do s = gr%s_i(), gr%s_o()

       pt%s = s

       k_i = gr%k_s_i(s)
       k_o = gr%k_s_o(s)

       cell_loop : do k = k_i, k_o-1

          pt%x = 0.5*(gr%pt(k)%x + gr%pt(k+1)%x)

          V_2 = this%coeff(I_V_2, pt)
          c_1 = this%coeff(I_C_1, pt)
          Gamma_1 = this%coeff(I_GAMMA_1, pt)

          I = I + SQRT(c_1*V_2/Gamma_1)*(gr%pt(k+1)%x - gr%pt(k)%x)

       end do cell_loop

    end do seg_loop
          
    Delta_p = 0.5_WP/I

    ! Finish

    return

  end function Delta_p

  !****

  function Delta_g (this, x_i, x_o, lambda)

    class(poly_model_t), intent(in) :: this
    real(WP), intent(in)            :: x_i
    real(WP), intent(in)            :: x_o
    real(WP), intent(in)            :: lambda
    real(WP)                        :: Delta_g

    type(grid_t)  :: gr
    real(WP)      :: I
    integer       :: s
    type(point_t) :: pt
    integer       :: k_i
    integer       :: k_o
    integer       :: k
    real(WP)      :: As
    real(WP)      :: c_1

    ! Evaluate the dimensionless g-mode inverse period separation,
    ! using a midpoint quadrature rule since the integrand typically
    ! diverges at the boundaries

    ! First, create the nested grid

    gr = grid_t(this%gr, x_i, x_o)

    ! Now evaluate the integrand segment by segment

    I = 0._WP

    seg_loop : do s = gr%s_i(), gr%s_o()

       pt%s = s

       k_i = gr%k_s_i(s)
       k_o = gr%k_s_o(s)

       cell_loop : do k = k_i, k_o-1

          pt%x = 0.5*(gr%pt(k)%x + gr%pt(k+1)%x)

          As = this%coeff(I_AS, pt)
          c_1 = this%coeff(I_C_1, pt)

          I = I + (SQRT(MAX(As/c_1, 0._WP))/pt%x)*(gr%pt(k+1)%x - gr%pt(k)%x)

       end do cell_loop

    end do seg_loop
          
    Delta_g = SQRT(lambda)/(2._WP*PI**2)*I

    ! Finish

    return

  end function Delta_g

  !****
  
  function grid (this) result (gr)

    class(poly_model_t), intent(in) :: this
    type(grid_t)                    :: gr

    ! Return the model grid

    gr = this%gr

    ! Finish

    return

  end function grid

end module gyre_poly_model
