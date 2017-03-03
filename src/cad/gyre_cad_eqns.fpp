! Module   : gyre_cad_eqns
! Purpose  : adiabatic differential equations (complex)
!
! Copyright 2017 Rich Townsend
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

module gyre_cad_eqns

  ! Uses

  use core_kinds

  use gyre_cad_trans
  use gyre_eqns
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: J_V_G = 1
  integer, parameter :: J_AS = 2
  integer, parameter :: J_U = 3
  integer, parameter :: J_C_1 = 4

  integer, parameter :: J_LAST = J_C_1

  ! Derived-type definitions

  type, extends (c_eqns_t) :: cad_eqns_t
     private
     class(model_t), pointer     :: ml => null()
     class(c_rot_t), allocatable :: rt
     type(cad_trans_t)           :: tr
     real(WP), allocatable       :: coeffs(:,:)
     real(WP), allocatable       :: x(:)
     real(WP)                    :: alpha_gr
     real(WP)                    :: alpha_om
   contains
     private
     procedure, public :: stencil
     procedure, public :: A
     procedure, public :: xA
  end type cad_eqns_t

  ! Interfaces

  interface cad_eqns_t
     module procedure cad_eqns_t_
  end interface cad_eqns_t

  ! Access specifiers

  private

  public :: cad_eqns_t

  ! Procedures

contains

  function cad_eqns_t_ (ml, pt_i, md_p, os_p) result (eq)

    class(model_t), pointer, intent(in) :: ml
    type(point_t), intent(in)           :: pt_i
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(cad_eqns_t)                    :: eq

    ! Construct the cad_eqns_t

    eq%ml => ml

    allocate(eq%rt, SOURCE=c_rot_t(ml, pt_i, md_p, os_p))

    eq%tr = ad_trans_t(ml, pt_i, md_p, os_p)

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

  end function cad_eqns_t_

  !****

  subroutine stencil (this, pt)

    class(cad_eqns_t), intent(inout) :: this
    type(point_t), intent(in)        :: pt(:)

    integer :: n_s
    integer :: i

    ! Calculate coefficients at the stencil points

    call check_model(this%ml, [I_V_2,I_AS,I_U,I_C_1,I_GAMMA_1])

    n_s = SIZE(pt)

    if (ALLOCATED(this%coeffs)) deallocate(this%coeffs)
    allocate(this%coeffs(n_s,J_LAST))

    do i = 1, n_s
       this%coeffs(i,J_V_G) = this%ml%coeff(I_V_2, pt(i))*pt(i)%x**2/this%ml%coeff(I_GAMMA_1, pt(i))
       this%coeffs(i,J_AS) = this%ml%coeff(I_AS, pt(i))
       this%coeffs(i,J_U) = this%ml%coeff(I_U, pt(i))
       this%coeffs(i,J_C_1) = this%ml%coeff(I_C_1, pt(i))
    end do

    this%x = pt%x

    ! Set up stencils for the rt and tr components

    call this%rt%stencil(pt)
    call this%tr%stencil(pt)

    ! Finish

    return

  end subroutine stencil

  !****

  function A (this, i, omega)

    class(cad_eqns_t), intent(in) :: this
    integer, intent(in)           :: i
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(i, omega)/this%x(i)

    ! Finish

    return

  end function A

  !****

  function xA (this, i, omega)

    class(cad_eqns_t), intent(in) :: this
    integer, intent(in)           :: i
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: xA(this%n_e,this%n_e)

    complex(WP) :: lambda
    complex(WP) :: l_i
    complex(WP) :: omega_c
    
    ! Evaluate the log(x)-space RHS matrix

    associate ( &
         V_g => this%coeffs(i,J_V_G), &
         As => this%coeffs(i,J_AS), &
         U => this%coeffs(i,J_U), &
         c_1 => this%coeffs(i,J_C_1), &
         alpha_gr => this%alpha_gr, &
         alpha_om => this%alpha_om)

      lambda = this%rt%lambda(i, omega)
      l_i = this%rt%l_i(omega)

      omega_c = this%rt%omega_c(i, omega)

      ! Set up the matrix

      xA(1,1) = V_g - 1._WP - l_i
      xA(1,2) = lambda/(c_1*alpha_om*omega_c**2) - V_g
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
      xA(4,2) = alpha_gr*(U*V_g)
      xA(4,3) = alpha_gr*(lambda)
      xA(4,4) = alpha_gr*(-U - l_i + 2._WP)

    end associate

    ! Apply the variables transformation

    call this%tr%trans_eqns(xA, i, omega)

    ! Finish

    return

  end function xA
  
end module gyre_ad_eqns
