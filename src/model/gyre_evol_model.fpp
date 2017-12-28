! Module   : gyre_evol_model
! Purpose  : stellar evolutionary model
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

module gyre_evol_model

  ! Uses

  use core_kinds
  
  use gyre_constants
  use gyre_grid
  use gyre_interp
  use gyre_model
  use gyre_model_par
  use gyre_point
  use gyre_util
  
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (model_t) :: evol_model_t
     private
     type(grid_t)                  :: gr
     type(r_interp_t), allocatable :: in(:,:)
     logical, allocatable          :: in_def(:)
     real(WP), public              :: M_star
     real(WP), public              :: R_star
     real(WP), public              :: L_star
     logical                       :: add_center
     logical                       :: repair_As
     integer                       :: s_i
     integer                       :: s_o
     character(:), allocatable     :: deriv_type
   contains
     private
     procedure, public :: define
     procedure, public :: coeff
     procedure, public :: dcoeff
     procedure, public :: M_r
     procedure, public :: P
     procedure, public :: rho
     procedure, public :: T
     procedure, public :: is_defined
     procedure, public :: is_vacuum
     procedure, public :: Delta_p
     procedure, public :: Delta_g
     procedure, public :: grid
  end type evol_model_t
 
  ! Interfaces

  interface evol_model_t
     module procedure evol_model_t_
  end interface evol_model_t

  ! Access specifiers

  private

  public :: evol_model_t

  ! Procedures

contains

  function evol_model_t_ (x, M_star, R_star, L_star, ml_p) result (ml)

    real(WP), intent(in)          :: x(:)
    real(WP), intent(in)          :: M_star
    real(WP), intent(in)          :: R_star
    real(WP), intent(in)          :: L_star
    type(model_par_t), intent(in) :: ml_p
    type(evol_model_t)            :: ml

    ! Construct the evol_model_t

    ! Create the grid

    if (ml_p%add_center) then

       if (x(1) /= 0._WP) then

          ml%gr = grid_t([0._WP,x])
          ml%add_center = .TRUE.

          if (check_log_level('INFO')) then
             write(OUTPUT_UNIT, 100) 'Added central point'
100          format(3X,A)
          endif

       else

          ml%gr = grid_t(x)
          ml%add_center = .FALSE.

          if (check_log_level('INFO')) then
             write(OUTPUT_UNIT, 100) 'No need to add central point'
          endif

       endif

    else

       ml%gr = grid_t(x)
       ml%add_center = .FALSE.

    endif

    ! Allocate arrays

    ml%s_i = ml%gr%s_i()
    ml%s_o = ml%gr%s_o()

    allocate(ml%in(I_LAST,ml%s_i:ml%s_o))

    allocate(ml%in_def(I_LAST))

    ml%in_def = .FALSE.

    ! Other initializations

    ml%M_star = M_star
    ml%R_star = R_star
    ml%L_star = L_star

    ml%repair_As = ml_p%repair_As
    ml%deriv_type = ml_p%deriv_type

    ! Finish

    return

  end function evol_model_t_

  !****

  subroutine define (this, i, coeff)

    class(evol_model_t), intent(inout) :: this
    integer, intent(in)                :: i
    real(WP), intent(in)               :: coeff(:)

    real(WP), allocatable :: coeff_(:)
    real(WP)              :: coeff_0
    integer               :: s

    $ASSERT_DEBUG(i >= 1 .AND. i <= I_LAST,Invalid index)

    ! Define the i'th coefficient

    ! If necessary, add a central point

    if (this%add_center) then

       select case (i)
       case (I_U)

          coeff_0 = 3._WP

       case (I_AS)

          coeff_0 = 0._WP

       case default

          ! Interpolate coeff at x=0 using parabolic fitting

          associate (x_1 => this%gr%pt(2)%x, &
                     x_2 => this%gr%pt(3)%x)
            coeff_0 = (x_2**2*coeff(1) - x_1**2*coeff(2))/(x_2**2 - x_1**2)
          end associate

       end select

       coeff_ = [coeff_0,coeff]

    else

       coeff_ = coeff

    endif

    $CHECK_BOUNDS(SIZE(coeff_),this%gr%n_k)

    ! If necessary, repair data at segment boundaries

    select case (i)
    case (I_AS)
       if (this%repair_As) call repair_coeff_(this%gr, coeff_)
    end select

    ! Set up per-segment interpolating splines
          
    seg_loop : do s = this%s_i, this%s_o

       associate (k_i => this%gr%k_s_i(s), k_o => this%gr%k_s_o(s))
         if (this%gr%pt(k_i)%x == 0._WP) then
            this%in(i, s) = r_interp_t(this%gr%pt(k_i:k_o)%x, coeff_(k_i:k_o), this%deriv_type, df_dx_a=0._WP)
         else
            this%in(i, s) = r_interp_t(this%gr%pt(k_i:k_o)%x, coeff_(k_i:k_o), this%deriv_type)
         endif
       end associate

    end do seg_loop

    this%in_def(i) = .TRUE.

    ! Finish

    return

  end subroutine define

  !****

  subroutine repair_coeff_ (gr, coeff)

    type(grid_t), intent(in) :: gr
    real(WP), intent(inout)  :: coeff(:)

    integer :: s_i
    integer :: s_o
    integer :: s
    integer :: k_i
    integer :: k_o
    
    ! Repair coefficient data at segment boundaries, via linear
    ! interpolation from the segment interior

    s_i = gr%s_i()
    s_o = gr%s_o()

    seg_loop : do s = s_i, s_o

       k_i = gr%k_s_i(s)
       k_o = gr%k_s_o(s)
         
         if (s > s_i .AND. k_i + 2 <= k_o) then
            coeff(k_i) = coeff(k_i+1) + (gr%pt(k_i)%x - gr%pt(k_i+1)%x)*(coeff(k_i+2) - coeff(k_i+1))/&
                 (gr%pt(k_i+2)%x - gr%pt(k_i+1)%x)
         endif
               
         if (s < s_o .AND. k_o - 2 >= k_i) then
            coeff(k_o) = coeff(k_o-1) + (gr%pt(k_o)%x - gr%pt(k_o-1)%x)*(coeff(k_o-1) - coeff(k_o-2))/ &
                 (gr%pt(k_o-1)%x - gr%pt(k_o-2)%x)
         endif

    end do seg_loop

    ! Finish

    return

  end subroutine repair_coeff_

  !****

  function coeff (this, i, pt)

    class(evol_model_t), intent(in) :: this
    integer, intent(in)             :: i
    type(point_t), intent(in)       :: pt
    real(WP)                        :: coeff

    $ASSERT_DEBUG(i >= 1 .AND. i <= I_LAST,Invalid index)
    $ASSERT_DEBUG(this%is_defined(i),Undefined coefficient)

    $ASSERT_DEBUG(pt%s >= this%s_i .AND. pt%s <= this%s_o,Invalid segment)

    ! Evaluate the i'th coefficient

    associate (s => pt%s, x => pt%x)
      coeff = this%in(i,s)%f(x)
    end associate

    ! Finish

    return

  end function coeff

  !****

  function dcoeff (this, i, pt)

    class(evol_model_t), intent(in) :: this
    integer, intent(in)             :: i
    type(point_t), intent(in)       :: pt
    real(WP)                        :: dcoeff

    $ASSERT_DEBUG(i >= 1 .AND. i <= I_LAST,Invalid index)
    $ASSERT_DEBUG(this%is_defined(i),Undefined coefficient)

    $ASSERT_DEBUG(pt%s >= this%s_i .AND. pt%s <= this%s_o,Invalid segment)

    ! Evaluate the logarithmic derivative of the i'th coefficient

    associate (s => pt%s, x => pt%x)
      if (x == 0._WP) then
         dcoeff = 0._WP
      else
         dcoeff = x*this%in(i,s)%df_dx(x)/this%in(i,s)%f(x)
      end if
    end associate

    ! Finish

    return

  end function dcoeff

  !****

  function M_r (this, pt)

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: M_r

    ! Evaluate the fractional mass coordinate

    M_r = this%M_star*(pt%x**3/this%coeff(I_C_1, pt))

    ! Finish

    return

  end function M_r
    
  !****

  function P (this, pt)

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: P

    ! Evaluate the total pressure

    P = (G_GRAVITY*this%M_star**2/(4._WP*PI*this%R_star**4))* &
        (this%coeff(I_U, pt)/(this%coeff(I_C_1, pt)**2*this%coeff(I_V_2, pt)))

    ! Finish

    return

  end function P
    
  !****

  function rho (this, pt)

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: rho

    ! Evaluate the density

    rho = (this%M_star/(4._WP*PI*this%R_star**3))*(this%coeff(I_U, pt)/this%coeff(I_C_1, pt))

    ! Finish

    return

  end function rho
    
  !****

  function T (this, pt)

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    real(WP)                        :: T

    ! Evaluate the temperature

    T = (3._WP*this%coeff(I_BETA_RAD, pt)*this%P(pt)/A_RADIATION)**0.25_WP

    ! Finish

    return

  end function T
    
  !****

  function is_defined (this, i)

    class(evol_model_t), intent(in) :: this
    integer, intent(in)             :: i
    logical                         :: is_defined

    $ASSERT_DEBUG(i >= 1 .AND. i <= I_LAST,Invalid index)

    ! Return the definition status of the i'th coefficient

    is_defined = this%in_def(i)

    ! Finish

    return

  end function is_defined

  !****

  function is_vacuum (this, pt)

    class(evol_model_t), intent(in) :: this
    type(point_t), intent(in)       :: pt
    logical                         :: is_vacuum

    $ASSERT_DEBUG(pt%s >= this%s_i .AND. pt%s <= this%s_o,Invalid segment)

    ! Return whether the point is a vacuum

    is_vacuum = this%coeff(I_U,pt) == 0._WP

    ! Finish

    return

  end function is_vacuum

  !****

  function Delta_p (this, x_i, x_o)

    class(evol_model_t), intent(in) :: this
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
    ! using a midpoint quadrature rule since the integrand can
    ! diverge at the surface

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

    class(evol_model_t), intent(in) :: this
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
    ! using a midpoint quadrature rule since the integrand can diverge
    ! at the boundaries

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

  function grid (this)

    class(evol_model_t), intent(in) :: this
    type(grid_t)                    :: grid

    ! Return the model grid

    grid = this%gr

    ! Finish

    return

  end function grid

end module gyre_evol_model
