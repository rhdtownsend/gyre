! Module   : gyre_nad_eqns
! Purpose  : nonadiabatic differential equations
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

module gyre_nad_eqns

  ! Uses

  use core_kinds

  use gyre_eqns
  use gyre_grid
  use gyre_linalg
  use gyre_mode_par
  use gyre_model
  use gyre_model_util
  use gyre_nad_vars
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: J_V = 1
  integer, parameter :: J_V_G = 2
  integer, parameter :: J_As = 3
  integer, parameter :: J_U = 4
  integer, parameter :: J_C_1 = 5
  integer, parameter :: J_DELTA = 6
  integer, parameter :: J_NABLA_AD = 7
  integer, parameter :: J_DNABLA_AD = 8
  integer, parameter :: J_NABLA = 9
  integer, parameter :: J_C_RAD = 10
  integer, parameter :: J_DC_RAD = 11
  integer, parameter :: J_C_THM = 12
  integer, parameter :: J_C_DIF = 13
  integer, parameter :: J_C_EPS_AD = 14
  integer, parameter :: J_C_EPS_S = 15
  integer, parameter :: J_KAP_AD = 16
  integer, parameter :: J_KAP_S = 17

  integer, parameter :: J_LAST = J_KAP_S

  ! Derived-type definitions

  type, extends (c_eqns_t) :: nad_eqns_t
     private
     class(model_t), pointer     :: ml => null()
     class(c_rot_t), allocatable :: rt
     type(nad_vars_t)            :: vr
     type(point_t), allocatable  :: pt(:)
     real(WP), allocatable       :: coeffs(:,:)
     real(WP)                    :: alpha_gr
     real(WP)                    :: alpha_hf
     real(WP)                    :: alpha_om
     logical                     :: narf_approx
   contains
     private
     procedure, public :: stencil
     procedure, public :: A
     procedure, public :: xA
  end type nad_eqns_t

  ! Interfaces

  interface nad_eqns_t
     module procedure nad_eqns_t_
  end interface nad_eqns_t

  ! Access specifiers

  private

  public :: nad_eqns_t

  ! Procedures

contains

  function nad_eqns_t_ (ml, gr, md_p, os_p) result (eq)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(nad_eqns_t)                    :: eq

    ! Construct the nad_eqns_t

    call check_model(ml, [ &
         I_V_2,I_AS,I_U,I_C_1,I_GAMMA_1,I_NABLA,I_NABLA_AD,I_DELTA,&
         I_C_RAD,I_C_THM,I_C_DIF,I_C_EPS_AD,I_C_EPS_S, &
         I_KAP_AD,I_KAP_S])

    eq%ml => ml

    allocate(eq%rt, SOURCE=c_rot_t(ml, gr, md_p, os_p))
    eq%vr = nad_vars_t(ml, gr, md_p, os_p)

    if (os_p%cowling_approx) then
       eq%alpha_gr = 0._WP
    else
       eq%alpha_gr = 1._WP
    endif

    if (os_p%narf_approx) then
       eq%alpha_hf = 0._WP
    else
       eq%alpha_hf = 1._WP
    endif

    select case (os_p%time_factor)
    case ('OSC')
       eq%alpha_om = 1._WP
    case ('EXP')
       eq%alpha_om = -1._WP
    case default
       $ABORT(Invalid time_factor)
    end select

    eq%n_e = 6

    ! Finish

    return

  end function nad_eqns_t_

  !****

  subroutine stencil (this, pt)

    class(nad_eqns_t), intent(inout) :: this
    type(point_t), intent(in)        :: pt(:)

    integer :: n_s
    integer :: i

    ! Calculate coefficients at the stencil points

    this%pt = pt

    n_s = SIZE(pt)

    if (ALLOCATED(this%coeffs)) deallocate(this%coeffs)
    allocate(this%coeffs(n_s, J_LAST))

    do i = 1, n_s
       this%coeffs(i,J_V) = this%ml%coeff(I_V_2, pt(i))*pt(i)%x**2
       this%coeffs(i,J_V_G) = this%coeffs(i,J_V)/this%ml%coeff(I_GAMMA_1, pt(i))
       this%coeffs(i,J_AS) = this%ml%coeff(I_AS, pt(i))
       this%coeffs(i,J_U) = this%ml%coeff(I_U, pt(i))
       this%coeffs(i,J_C_1) = this%ml%coeff(I_C_1, pt(i))
       this%coeffs(i,J_NABLA_AD) = this%ml%coeff(I_NABLA_AD, pt(i))
       this%coeffs(i,J_DNABLA_AD) = this%ml%dcoeff(I_NABLA_AD, pt(i))
       this%coeffs(i,J_NABLA) = this%ml%coeff(I_NABLA, pt(i))
       this%coeffs(i,J_DELTA) = this%ml%coeff(I_DELTA, pt(i))
       this%coeffs(i,J_C_RAD) = this%ml%coeff(I_C_RAD, pt(i))
       this%coeffs(i,J_DC_RAD) = this%ml%dcoeff(I_C_RAD, pt(i))
       this%coeffs(i,J_C_THM) = this%ml%coeff(I_C_THM, pt(i))
       this%coeffs(i,J_C_DIF) = this%ml%coeff(I_C_DIF, pt(i))
       this%coeffs(i,J_C_EPS_AD) = this%ml%coeff(I_C_EPS_AD, pt(i))
       this%coeffs(i,J_C_EPS_S) = this%ml%coeff(I_C_EPS_S, pt(i))
       this%coeffs(i,J_KAP_AD) = this%ml%coeff(I_KAP_AD, pt(i))
       this%coeffs(i,J_KAP_S) = this%ml%coeff(I_KAP_S, pt(i))
    end do

    ! Finish

    return

  end subroutine stencil

  !****

  function A (this, i, omega)

    class(nad_eqns_t), intent(in) :: this
    integer, intent(in)           :: i
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(i, omega)/this%pt(i)%x

    ! Finish

    return

  end function A

  !****

  function xA (this, i, omega)

    class(nad_eqns_t), intent(in) :: this
    integer, intent(in)           :: i
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: xA(this%n_e,this%n_e)

    complex(WP) :: lambda
    complex(WP) :: l_i
    complex(WP) :: omega_c
         
    ! Evaluate the log(x)-space RHS matrix

    associate ( &
         V => this%coeffs(i,J_V), &
         V_g => this%coeffs(i,J_V_G), &
         As => this%coeffs(i,J_AS), &
         U => this%coeffs(i,J_U), &
         c_1 => this%coeffs(i,J_C_1), &
         nabla_ad => this%coeffs(i,J_NABLA_AD), &
         dnabla_ad => this%coeffs(i,J_DNABLA_AD), &
         nabla => this%coeffs(i,J_NABLA), &
         delta => this%coeffs(i,J_DELTA), &
         c_rad => this%coeffs(i,J_C_RAD), &
         dc_rad => this%coeffs(i,J_DC_RAD), &
         c_thm => this%coeffs(i,J_C_THM), &
         c_dif => this%coeffs(i,J_C_DIF), &
         c_eps_ad => this%coeffs(i,J_C_EPS_AD), &
         c_eps_S => this%coeffs(i,J_C_EPS_S), &
         kap_ad => this%coeffs(i,J_KAP_AD), &
         kap_S => this%coeffs(i,J_KAP_S), &
         alpha_gr => this%alpha_gr, &
         alpha_hf => this%alpha_hf, &
         alpha_om => this%alpha_om, &
         pt => this%pt(i))
         
      lambda = this%rt%lambda(pt, omega)
      l_i = this%rt%l_i(omega)
    
      omega_c = this%rt%omega_c(pt, omega)

      ! Set up the matrix

      xA(1,1) = V_g - 1._WP - l_i
      xA(1,2) = lambda/(c_1*alpha_om*omega_c**2) - V_g
      xA(1,3) = alpha_gr*(lambda/(c_1*alpha_om*omega_c**2))
      xA(1,4) = alpha_gr*(0._WP)
      xA(1,5) = delta
      xA(1,6) = 0._WP

      xA(2,1) = c_1*alpha_om*omega_c**2 - As
      xA(2,2) = As - U + 3._WP - l_i
      xA(2,3) = alpha_gr*(0._WP)
      xA(2,4) = alpha_gr*(-1._WP)
      xA(2,5) = delta
      xA(2,6) = 0._WP

      xA(3,1) = alpha_gr*(0._WP)
      xA(3,2) = alpha_gr*(0._WP)
      xA(3,3) = alpha_gr*(3._WP - U - l_i)
      xA(3,4) = alpha_gr*(1._WP)
      xA(3,5) = alpha_gr*(0._WP)
      xA(3,6) = alpha_gr*(0._WP)

      xA(4,1) = alpha_gr*(U*As)
      xA(4,2) = alpha_gr*(U*V_g)
      xA(4,3) = alpha_gr*(lambda)
      xA(4,4) = alpha_gr*(-U - l_i + 2._WP)
      xA(4,5) = alpha_gr*(-U*delta)
      xA(4,6) = alpha_gr*(0._WP)

      xA(5,1) = V*(nabla_ad*(U - c_1*alpha_om*omega_c**2) - 4._WP*(nabla_ad - nabla) + c_dif + nabla_ad*dnabla_ad)
      xA(5,2) = V*(lambda/(c_1*alpha_om*omega_c**2)*(nabla_ad - nabla) - (c_dif + nabla_ad*dnabla_ad))
      xA(5,3) = alpha_gr*(V*lambda/(c_1*alpha_om*omega_c**2)*(nabla_ad - nabla))
      xA(5,4) = alpha_gr*(V*nabla_ad)
      xA(5,5) = V*nabla*(4._WP - kap_S) - (l_i - 2._WP)
      xA(5,6) = -V*nabla/c_rad

      xA(6,1) = alpha_hf*lambda*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad
      xA(6,2) = V*c_eps_ad - lambda*c_rad*(alpha_hf*nabla_ad/nabla - (3._WP + dc_rad)/(c_1*alpha_om*omega_c**2))
      xA(6,3) = alpha_gr*(lambda*c_rad*(3._WP + dc_rad)/(c_1*alpha_om*omega_c**2))
      xA(6,4) = alpha_gr*(0._WP)
      if (pt%x > 0._WP) then
         xA(6,5) = c_eps_S - alpha_hf*lambda*c_rad/(nabla*V) + (0._WP,1._WP)*SQRT(CMPLX(alpha_om, KIND=WP))*omega_c*c_thm
      else
         xA(6,5) = -alpha_hf*HUGE(0._WP)
      endif
      xA(6,6) = -1._WP - l_i

      ! Apply the variables transformation

      xA = MATMUL(this%vr%G(pt, omega), MATMUL(xA, this%vr%H(pt, omega)) - &
                                                   this%vr%dH(pt, omega))

    end associate

    ! Finish

    return

  end function xA

end module gyre_nad_eqns
