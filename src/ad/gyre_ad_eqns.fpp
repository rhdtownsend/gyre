! Module   : gyre_ad_eqns
! Purpose  : adiabatic differential equations
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
  integer, parameter :: J_OMEGA_ROT = 6
  integer, parameter :: J_OMEGA_ROT_I = 7

  integer, parameter :: J_LAST = J_OMEGA_ROT_I

  ! Derived-type definitions

  type, extends (r_eqns_t) :: ad_eqns_t
     private
     type(context_t), pointer :: cx => null()
     type(ad_trans_t)         :: tr
     real(WP), allocatable    :: coeff(:,:)
     real(WP), allocatable    :: x(:)
     real(WP)                 :: alpha_gr
     real(WP)                 :: alpha_om
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

    if (os_p%cowling_approx) then
       eq%alpha_gr = 0._WP
    else
       eq%alpha_gr = 1._WP
    endif

    select case (os_p%time_factor)
    case ('OSC')
       eq%alpha_om = 1._WP
    case ('EXP')
       eq%alpha_om = -1._WP
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

    integer :: n_s
    integer :: i

    ! Calculate coefficients at the stencil points

    associate (ml => this%cx%ml)

      call check_model(ml, [I_V_2,I_AS,I_U,I_C_1,I_GAMMA_1,I_OMEGA_ROT])

      n_s = SIZE(pt)

      if (ALLOCATED(this%coeff)) deallocate(this%coeff)
      allocate(this%coeff(n_s,J_LAST))

      do i = 1, n_s
         this%coeff(i,J_V) = ml%coeff(I_V_2, pt(i))*pt(i)%x**2
         this%coeff(i,J_AS) = ml%coeff(I_AS, pt(i))
         this%coeff(i,J_U) = ml%coeff(I_U, pt(i))
         this%coeff(i,J_C_1) = ml%coeff(I_C_1, pt(i))
         this%coeff(i,J_GAMMA_1) = ml%coeff(I_GAMMA_1, pt(i))
         this%coeff(i,J_OMEGA_ROT) = ml%coeff(I_OMEGA_ROT, pt(i))
      end do

      this%coeff(:,J_OMEGA_ROT_I) = ml%coeff(I_OMEGA_ROT, this%cx%pt_i)

      this%x = pt%x

    end associate

    ! Set up stencil for the tr component

    call this%tr%stencil(pt)

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

    A = this%xA(i, st)/this%x(i)

    ! Finish

    return

  end function A

  !****

  function xA (this, i, st)

    class(ad_eqns_t), intent(in) :: this
    integer, intent(in)          :: i
    class(r_state_t), intent(in) :: st
    real(WP)                     :: xA(this%n_e,this%n_e)

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
         Omega_rot => this%coeff(i,J_OMEGA_ROT), &
         Omega_rot_i => this%coeff(i,J_OMEGA_ROT_I), &
         alpha_gr => this%alpha_gr, &
         alpha_om => this%alpha_om)

      omega_c = this%cx%omega_c(Omega_rot, st)

      lambda = this%cx%lambda(Omega_rot, st)
      l_i = this%cx%l_e(Omega_rot_i, st)

      ! Set up the matrix

      xA(1,1) = V/Gamma_1 - 1._WP - l_i
      xA(1,2) = lambda/(c_1*alpha_om*omega_c**2) - V/Gamma_1
      xA(1,3) = alpha_gr*(lambda/(c_1*alpha_om*omega_c**2))
      xA(1,4) = alpha_gr*(0._WP)

      xA(2,1) = c_1*alpha_om*omega_c**2 - As
      xA(2,2) = As - U + 3._WP - l_i
      xA(2,3) = alpha_gr*(0._WP)
      xA(2,4) = alpha_gr*(-1._WP)

      xA(3,1) = alpha_gr*(0._WP)
      xA(3,2) = alpha_gr*(0._WP)
      xA(3,3) = alpha_gr*(3._WP - U - l_i)
      xA(3,4) = alpha_gr*(1._WP)

      xA(4,1) = alpha_gr*(U*As)
      xA(4,2) = alpha_gr*(U*V/Gamma_1)
      xA(4,3) = alpha_gr*(lambda)
      xA(4,4) = alpha_gr*(-U - l_i + 2._WP)

    end associate

    ! Apply the variables transformation

    call this%tr%trans_eqns(xA, i, st)

    ! Finish

    return

  end function xA
  
end module gyre_ad_eqns
