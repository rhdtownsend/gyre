! Module   : gyre_nad_eqns
! Purpose  : nonadiabatic differential equations
!
! Copyright 2013-2020 Rich Townsend & The GYRE Team
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

  use gyre_context
  use gyre_eqns
  use gyre_linalg
  use gyre_math
  use gyre_model
  use gyre_mode_par
  use gyre_model_util
  use gyre_nad_trans
  use gyre_osc_par
  use gyre_point
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: P1_CONV_SCHEME = 1
  integer, parameter :: P4_CONV_SCHEME = 4

  integer, parameter :: J_V = 1
  integer, parameter :: J_As = 2
  integer, parameter :: J_U = 3
  integer, parameter :: J_C_1 = 4
  integer, parameter :: J_GAMMA_1 = 5
  integer, parameter :: J_DELTA = 6
  integer, parameter :: J_NABLA_AD = 7
  integer, parameter :: J_DNABLA_AD = 8
  integer, parameter :: J_NABLA = 9
  integer, parameter :: J_C_LUM = 10
  integer, parameter :: J_DC_LUM = 11
  integer, parameter :: J_C_RAD = 12
  integer, parameter :: J_DC_RAD = 13
  integer, parameter :: J_C_THN = 14
  integer, parameter :: J_DC_THN = 15
  integer, parameter :: J_C_THK = 16
  integer, parameter :: J_C_EPS = 17
  integer, parameter :: J_KAP_RHO = 18
  integer, parameter :: J_KAP_T = 19

  integer, parameter :: J_LAST = J_KAP_T

  ! Derived-type definitions

  type, extends (c_eqns_t) :: nad_eqns_t
     private
     type(context_t), pointer   :: cx => null()
     type(point_t), allocatable :: pt(:)
     type(nad_trans_t)          :: tr
     real(WP), allocatable      :: coeff(:,:)
     real(WP)                   :: x_atm
     real(WP)                   :: gamma_gr
     real(WP)                   :: gamma_th
     real(WP)                   :: gamma_hf
     real(WP)                   :: gamma_rh
     real(WP)                   :: gamma_om
     real(WP)                   :: alpha_pi
     real(WP)                   :: alpha_gamma
     integer                    :: conv_scheme
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

  function nad_eqns_t_ (cx, md_p, os_p) result (eq)

    type(context_t), pointer, intent(in) :: cx
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(nad_eqns_t)                     :: eq
    type(point_t)                        :: pt_o

    ! Construct the nad_eqns_t

    eq%cx => cx

    eq%tr = nad_trans_t(cx, md_p, os_p)

    eq%gamma_gr = os_p%gamma_gr
    eq%gamma_th = os_p%gamma_th
    eq%gamma_hf = os_p%gamma_hf
       
    if (os_p%eddington_approx) then
       eq%gamma_rh = 1._WP
    else
       eq%gamma_rh = 0._WP
    endif

    select case (os_p%time_factor)
    case ('OSC')
       eq%gamma_om = 1._WP
    case ('EXP')
       eq%gamma_om = -1._WP
    case default
       $ABORT(Invalid time_factor)
    end select

    select case (os_p%conv_scheme)
    case ('FROZEN_PESNELL_1')
       eq%conv_scheme = P1_CONV_SCHEME
    case ('FROZEN_PESNELL_4')
       eq%conv_scheme = P4_CONV_SCHEME
    case default
       $ABORT(Invalid conv_scheme)
    end select

    eq% x_atm = os_p%x_atm
    if (eq% x_atm < 0._WP) then
       pt_o = cx% point_o()
       eq% x_atm = pt_o% x
    end if

    select case (os_p%isolation)
    case ('GAMMA')
       eq%alpha_gamma = 0._WP
       eq%alpha_pi = 1._WP
    case ('PI')
       eq%alpha_gamma = 1._WP
       eq%alpha_pi = 0._WP
    case ('NONE')
       eq%alpha_gamma = 1._WP
       eq%alpha_pi = 1._WP
    case default
       $ABORT(Invalid isolation condition)
    end select

    eq%n_e = 6

    ! Finish

    return

  end function nad_eqns_t_

  !****

  subroutine stencil (this, pt)

    class(nad_eqns_t), intent(inout) :: this
    type(point_t), intent(in)        :: pt(:)

    class(model_t), pointer :: ml
    integer                 :: n_s
    integer                 :: i

    ! Calculate coefficients at the stencil points

    ml => this%cx%model()

    call check_model(ml, [ &
         I_V_2,I_AS,I_U,I_C_1,I_GAMMA_1,I_NABLA,I_NABLA_AD,I_DELTA, &
         I_C_LUM,I_C_RAD,I_C_THN,I_C_THK,I_C_EPS, &
         I_KAP_RHO,I_KAP_T])

    n_s = SIZE(pt)

    if (ALLOCATED(this%coeff)) deallocate(this%coeff)
    allocate(this%coeff(n_s, J_LAST))

    do i = 1, n_s

       $ASSERT(.NOT. ml%is_vacuum(pt(i)),Attempt to stencil at vacuum point)

       this%coeff(i,J_V) = ml%coeff(I_V_2, pt(i))*pt(i)%x**2
       this%coeff(i,J_AS) = ml%coeff(I_AS, pt(i))
       this%coeff(i,J_U) = ml%coeff(I_U, pt(i))
       this%coeff(i,J_C_1) = ml%coeff(I_C_1, pt(i))
       this%coeff(i,J_GAMMA_1) = ml%coeff(I_GAMMA_1, pt(i))
       this%coeff(i,J_NABLA_AD) = ml%coeff(I_NABLA_AD, pt(i))
       this%coeff(i,J_DNABLA_AD) = ml%dcoeff(I_NABLA_AD, pt(i))
       this%coeff(i,J_NABLA) = ml%coeff(I_NABLA, pt(i))
       this%coeff(i,J_DELTA) = ml%coeff(I_DELTA, pt(i))
       this%coeff(i,J_C_LUM) = ml%coeff(I_C_LUM, pt(i))
       this%coeff(i,J_DC_LUM) = ml%dcoeff(I_C_LUM, pt(i))
       this%coeff(i,J_C_RAD) = ml%coeff(I_C_RAD, pt(i))
       this%coeff(i,J_DC_RAD) = ml%dcoeff(I_C_RAD, pt(i))
       this%coeff(i,J_C_THN) = ml%coeff(I_C_THN, pt(i))
       this%coeff(i,J_DC_THN) = ml%dcoeff(I_C_THN, pt(i))
       this%coeff(i,J_C_THK) = ml%coeff(I_C_THK, pt(i))
       this%coeff(i,J_C_EPS) = ml%coeff(I_C_EPS, pt(i))
       this%coeff(i,J_KAP_RHO) = ml%coeff(I_KAP_RHO, pt(i))
       this%coeff(i,J_KAP_T) = ml%coeff(I_KAP_T, pt(i))

    end do

    ! Set up stencil for the tr component

    call this%tr%stencil(pt)

    ! Store the stencil points for on-the-fly evaluations

    this%pt = pt

    ! Finish

    return

  end subroutine stencil

  !****

  function A (this, i, st)

    class(nad_eqns_t), intent(in) :: this
    integer, intent(in)           :: i
    class(c_state_t), intent(in)  :: st
    complex(WP)                   :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(i, st)/this%pt(i)%x

    ! Finish

    return

  end function A

  !****

  function xA (this, i, st)

    class(nad_eqns_t), intent(in) :: this
    integer, intent(in)           :: i
    class(c_state_t), intent(in)  :: st
    complex(WP)                   :: xA(this%n_e,this%n_e)

    real(WP)    :: Omega_rot
    real(WP)    :: Omega_rot_i
    complex(WP) :: omega_c
    complex(WP) :: i_omega_c
    complex(WP) :: lambda
    complex(WP) :: l_i
    complex(WP) :: f_rh
    complex(WP) :: df_rh
    complex(WP) :: conv_term
    complex(WP) :: eps_rho
    complex(WP) :: eps_T
    complex(WP) :: c_eps_ad
    complex(WP) :: c_eps_S
    real(WP)    :: kap_ad
    real(WP)    :: kap_S
    real(WP)    :: c_dif
         
    ! Evaluate the log(x)-space RHS matrix

    associate ( &
         V => this%coeff(i,J_V), &
         As => this%coeff(i,J_AS), &
         U => this%coeff(i,J_U), &
         c_1 => this%coeff(i,J_C_1), &
         Gamma_1 => this%coeff(i,J_GAMMA_1), &
         nabla_ad => this%coeff(i,J_NABLA_AD), &
         dnabla_ad => this%coeff(i,J_DNABLA_AD), &
         nabla => this%coeff(i,J_NABLA), &
         delta => this%coeff(i,J_DELTA), &
         c_lum => this%coeff(i,J_C_LUM), &
         dc_lum => this%coeff(i,J_DC_LUM), &
         c_rad => this%coeff(i,J_C_RAD), &
         dc_rad => this%coeff(i,J_DC_RAD), &
         c_thn => this%coeff(i,J_C_THN), &
         dc_thn => this%coeff(i,J_DC_THN), &
         c_thk => this%coeff(i,J_C_THK), &
         c_eps => this%coeff(i,J_C_EPS), &
         kap_rho => this%coeff(i,J_KAP_RHO), &
         kap_T => this%coeff(i,J_KAP_T), &
         pt => this%pt(i), &
         pt_i => this%cx%point_i(), &
         x => this%pt(i)%x, &
         x_atm => this%x_atm, &
         gamma_gr => this%gamma_gr, &
         gamma_th => this%gamma_th, &
         gamma_hf => this%gamma_hf, &
         gamma_rh => this%gamma_rh, &
         gamma_om => this%gamma_om, &
         alpha_pi => this%alpha_pi, &
         alpha_gamma => this%alpha_gamma)

      Omega_rot = this%cx%Omega_rot(pt)
      Omega_rot_i = this%cx%Omega_rot(pt_i)

      omega_c = this%cx%omega_c(Omega_rot, st)
      omega_c = this%cx%omega_c(Omega_rot, st)
      i_omega_c = (0._WP,1._WP)*sqrt(CMPLX(gamma_om, KIND=WP))*omega_c

      lambda = this%cx%lambda(Omega_rot, st)
      l_i = this%cx%l_e(Omega_rot_i, st)
    
      f_rh = 1._WP - 0.25_WP*gamma_rh*i_omega_c*c_thn
      df_rh = -0.25_WP*gamma_rh*i_omega_c*c_thn*dc_thn/f_rh

      select case (this%conv_scheme)
      case (P1_CONV_SCHEME)
         conv_term = lambda*c_rad*(3._WP + dc_rad)/(c_1*gamma_om*omega_c**2)
      case (P4_CONV_SCHEME)
         conv_term = lambda*(c_lum*(3._WP + dc_lum) - (c_lum - c_rad))/(c_1*gamma_om*omega_c**2)
      case default
         $ABORT(Invalid conv_scheme)
      end select

      eps_rho = this%cx%eps_rho(st, pt)
      eps_T = this%cx%eps_T(st, pt)

      c_eps_ad = c_eps*(nabla_ad*eps_T + eps_rho/Gamma_1)
      c_eps_S = c_eps*(eps_T - delta*eps_rho)

      kap_ad = nabla_ad*kap_T + kap_rho/Gamma_1
      kap_S = kap_T - delta*kap_rho
      
      c_dif = (kap_ad-4._WP*nabla_ad)*V*nabla + nabla_ad*(dnabla_ad + V)

      ! Set up the matrix

      xA(1,1) = V/Gamma_1 - 1._WP - l_i
      xA(1,2) = lambda/(c_1*gamma_om*omega_c**2) - V/Gamma_1 * alpha_gamma
      xA(1,3) = gamma_gr*(lambda/(c_1*gamma_om*omega_c**2))
      xA(1,4) = gamma_gr*(0._WP)
      xA(1,5) = delta
      xA(1,6) = 0._WP

      xA(2,1) = c_1*gamma_om*omega_c**2 - As * MERGE(MERGE(alpha_pi, alpha_gamma, x<x_atm), 1._WP, As > 0)
      xA(2,2) = As - U + 3._WP - l_i
      xA(2,3) = gamma_gr*(0._WP)
      xA(2,4) = gamma_gr*(-1._WP)
      xA(2,5) = delta
      xA(2,6) = 0._WP

      xA(3,1) = gamma_gr*(0._WP)
      xA(3,2) = gamma_gr*(0._WP)
      xA(3,3) = gamma_gr*(3._WP - U - l_i)
      xA(3,4) = gamma_gr*(1._WP)
      xA(3,5) = gamma_gr*(0._WP)
      xA(3,6) = gamma_gr*(0._WP)

      xA(4,1) = gamma_gr*(U*As)
      xA(4,2) = gamma_gr*(U*V/Gamma_1)
      xA(4,3) = gamma_gr*(lambda)
      xA(4,4) = gamma_gr*(-U - l_i + 2._WP)
      xA(4,5) = gamma_gr*(-U*delta)
      xA(4,6) = gamma_gr*(0._WP)

      xA(5,1) = V*(nabla_ad*(U - c_1*gamma_om*omega_c**2) - 4._WP*(nabla_ad - nabla) + c_dif)/f_rh
      xA(5,2) = V*(lambda/(c_1*gamma_om*omega_c**2)*(nabla_ad - nabla) - c_dif)/f_rh
      xA(5,3) = gamma_gr*(V*lambda/(c_1*gamma_om*omega_c**2)*(nabla_ad - nabla))/f_rh
      xA(5,4) = gamma_gr*(V*nabla_ad)/f_rh
      xA(5,5) = V*nabla*(4._WP*f_rh - kap_S)/f_rh - df_rh - (l_i - 2._WP)
      xA(5,6) = -V*nabla/(c_rad*f_rh)

      xA(6,1) = gamma_hf*lambda*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad
      xA(6,2) = V*c_eps_ad - lambda*c_rad*gamma_hf*nabla_ad/nabla + conv_term
      xA(6,3) = gamma_gr*conv_term
      xA(6,4) = gamma_gr*(0._WP)
      if (x > 0._WP) then
         xA(6,5) = c_eps_S - gamma_hf*lambda*c_rad/(nabla*V) + gamma_th*i_omega_c*c_thk
      else
         xA(6,5) = -gamma_hf*HUGE(0._WP)
      endif
      xA(6,6) = -1._WP - l_i

    end associate

    ! Apply the variables transformation

    call this%tr%trans_eqns(xA, i, st)

    ! Finish

    return

  end function xA

end module gyre_nad_eqns
