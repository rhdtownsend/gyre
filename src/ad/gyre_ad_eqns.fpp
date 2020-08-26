! Module   : gyre_ad_eqns
! Purpose  : adiabatic differential equations
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

module gyre_ad_eqns

  ! Uses

  use core_kinds

  use gyre_ad_trans
  use gyre_context
  use gyre_eqns
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: J_V = 1
  integer, parameter :: J_AS = 2
  integer, parameter :: J_U = 3
  integer, parameter :: J_C_1 = 4
  integer, parameter :: J_GAMMA_1 = 5

  integer, parameter :: J_LAST = J_GAMMA_1

  ! Derived-type definitions

  type, extends (r_eqns_t) :: ad_eqns_t
     private
     type(context_t), pointer   :: cx => null()
     type(point_t), allocatable :: pt(:)
     type(ad_trans_t)           :: tr
     real(WP), allocatable      :: coeff(:,:)
     real(WP)                   :: x_atm
     real(WP)                   :: gamma_gr
     real(WP)                   :: alpha_pi
     real(WP)                   :: alpha_gamma
     real(WP)                   :: gamma_om
   contains
     private
     procedure, public :: stencil
     procedure, public :: A
     procedure, public :: xA
  end type ad_eqns_t

  ! Interfaces

  interface ad_eqns_t
     module procedure ad_eqns_t_
  end interface ad_eqns_t

  ! Access specifiers

  private

  public :: ad_eqns_t

  ! Procedures

contains

  function ad_eqns_t_ (cx, md_p, os_p) result (eq)

    type(context_t), pointer, intent(in) :: cx
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(ad_eqns_t)                      :: eq

    ! Construct the ad_eqns_t

    eq%cx => cx

    eq%tr = ad_trans_t(cx, md_p, os_p)

    eq%gamma_gr = os_p%gamma_gr

    eq%x_atm = os_p%x_atm
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

    select case (os_p%time_factor)
    case ('OSC')
       eq%gamma_om = 1._WP
    case ('EXP')
       eq%gamma_om = -1._WP
    case default
       $ABORT(Invalid time_factor)
    end select

    eq%n_e = 4

    ! Finish

    return

  end function ad_eqns_t_

  !****

  subroutine stencil (this, pt)

    class(ad_eqns_t), intent(inout) :: this
    type(point_t), intent(in)       :: pt(:)

    class(model_t), pointer :: ml
    integer                 :: n_s
    integer                 :: i

    ! Calculate coefficients at the stencil points

    ml => this%cx%model()

    call check_model(ml, [I_V_2,I_AS,I_U,I_C_1,I_GAMMA_1])

    n_s = SIZE(pt)

    if (ALLOCATED(this%coeff)) deallocate(this%coeff)
    allocate(this%coeff(n_s,J_LAST))

    do i = 1, n_s

       $ASSERT(.NOT. ml%is_vacuum(pt(i)),Attempt to stencil at vacuum point)
       
       this%coeff(i,J_V) = ml%coeff(I_V_2, pt(i))*pt(i)%x**2
       this%coeff(i,J_AS) = ml%coeff(I_AS, pt(i))
       this%coeff(i,J_U) = ml%coeff(I_U, pt(i))
       this%coeff(i,J_C_1) = ml%coeff(I_C_1, pt(i))
       this%coeff(i,J_GAMMA_1) = ml%coeff(I_GAMMA_1, pt(i))

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

    class(ad_eqns_t), intent(in) :: this
    integer, intent(in)          :: i
    class(r_state_t), intent(in) :: st
    real(WP)                     :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(i, st)/this%pt(i)%x

    ! Finish

    return

  end function A

  !****

  function xA (this, i, st)

    class(ad_eqns_t), intent(in) :: this
    integer, intent(in)          :: i
    class(r_state_t), intent(in) :: st
    real(WP)                     :: xA(this%n_e,this%n_e)

    real(WP) :: Omega_rot
    real(WP) :: Omega_rot_i
    real(WP) :: omega_c
    real(WP) :: lambda
    real(WP) :: l_i
    
    ! Evaluate the log(x)-space RHS matrix

    associate ( &
         V => this%coeff(i,J_V), &
         As => this%coeff(i,J_AS), &
         U => this%coeff(i,J_U), &
         c_1 => this%coeff(i,J_C_1), &
         Gamma_1 => this%coeff(i,J_GAMMA_1), &
         pt => this%pt(i), &
         pt_i => this%cx%point_i(), &
         alpha_pi => this%alpha_pi, &
         alpha_gamma => this%alpha_gamma, &
         x => this%pt(i)%x, &
         x_atm => this%x_atm, &
         gamma_gr => this%gamma_gr, &
         gamma_om => this%gamma_om)

      Omega_rot = this%cx%Omega_rot(pt)
      Omega_rot_i = this%cx%Omega_rot(pt_i)

      omega_c = this%cx%omega_c(Omega_rot, st)

      lambda = this%cx%lambda(Omega_rot, st)
      l_i = this%cx%l_e(Omega_rot_i, st)

      ! Set up the matrix

      xA(1,1) = V/Gamma_1 - 1._WP - l_i
      xA(1,2) = lambda/(c_1*gamma_om*omega_c**2) - V/Gamma_1 * alpha_gamma
      xA(1,3) = gamma_gr*(lambda/(c_1*gamma_om*omega_c**2))
      xA(1,4) = gamma_gr*(0._WP)

      xA(2,1) = c_1*gamma_om*omega_c**2 - As * MERGE(MERGE(alpha_pi, alpha_gamma, x<x_atm), 1._WP, As > 0)
      xA(2,2) = As - U + 3._WP - l_i
      xA(2,3) = gamma_gr*(0._WP)
      xA(2,4) = gamma_gr*(-1._WP)

      xA(3,1) = gamma_gr*(0._WP)
      xA(3,2) = gamma_gr*(0._WP)
      xA(3,3) = gamma_gr*(3._WP - U - l_i)
      xA(3,4) = gamma_gr*(1._WP)

      xA(4,1) = gamma_gr*(U*As)
      xA(4,2) = gamma_gr*(U*V/Gamma_1)
      xA(4,3) = gamma_gr*(lambda)
      xA(4,4) = gamma_gr*(-U - l_i + 2._WP)

    end associate

    ! Apply the variables transformation

    call this%tr%trans_eqns(xA, i, st)

    ! Finish

    return

  end function xA
  
end module gyre_ad_eqns
