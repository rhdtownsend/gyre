! Module   : gyre_rad_vars
! Purpose  : adiabatic radial variables transformations
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_rad_vars

  ! Uses

  use core_kinds

  use gyre_linalg
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: DZIEM_SET = 1
  integer, parameter :: JCD_SET = 2
  integer, parameter :: MIX_SET = 3
  integer, parameter :: LAGP_SET = 4

  ! Derived-type definitions

  type :: rad_vars_t
     private
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     integer                     :: set
   contains
     private
     procedure, public :: G
     procedure         :: G_jcd_
     procedure         :: G_lagp_
     procedure, public :: H
     procedure         :: H_jcd_
     procedure         :: H_lagp_
     procedure, public :: dH
     procedure         :: dH_jcd_
     procedure         :: dH_lagp_
  end type rad_vars_t

  ! Interfaces

  interface rad_vars_t
     module procedure rad_vars_t_
  end interface rad_vars_t

  ! Access specifiers

  private

  public :: rad_vars_t

  ! Procedures

contains

  function rad_vars_t_ (ml, md_p, os_p) result (vr)

    class(model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(rad_vars_t)                    :: vr

    ! Construct the rad_vars_t

    vr%ml => ml

    allocate(vr%rt, SOURCE=r_rot_t(ml, md_p, os_p))

    select case (os_p%variables_set)
    case ('DZIEM')
       vr%set = DZIEM_SET
    case ('JCD')
       vr%set = JCD_SET
    case ('MIX')
       vr%set = MIX_SET
    case ('LAGP')
       vr%set = LAGP_SET
    case default
       $ABORT(Invalid variables_set)
    end select

    ! Finish

    return

  end function rad_vars_t_

  !****

  function G (this, s, x, omega)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: s
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: G(2,2)

    ! Evaluate the transformation matrix to convert variables from
    ! the canonical form

    select case (this%set)
    case (DZIEM_SET)
       G = identity_matrix(2)
    case (JCD_SET)
       G = this%G_jcd_(s, x, omega)
    case (MIX_SET)
       G = identity_matrix(2)
    case (LAGP_SET)
       G = this%G_lagp_(s, x, omega)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return
    
  end function G
  
  !****

  function G_jcd_ (this, s, x, omega) result (G)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: s
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: G(2,2)

    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! from the canonical form

    ! Calculate coefficients

    U = this%ml%U(s, x)
    c_1 = this%ml%c_1(s, x)

    omega_c = this%rt%omega_c(s, x, omega)

    ! Set up the matrix

    G(1,1) = 1._WP
    G(1,2) = 0._WP

    G(2,1) = 0._WP
    G(2,2) = 1._WP/(c_1*omega_c**2)

    ! Finish

    return

  end function G_jcd_

  !****

  function G_lagp_ (this, s, x, omega) result (G)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: s
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: G(2,2)

    real(WP) :: V_2

    ! Evaluate the transformation matrix to convert LAGP variables
    ! from the canonical form

    ! Calculate coefficients

    V_2 = this%ml%V_2(s, x)

    ! Set up the matrix

    G(1,1) = 1._WP
    G(1,2) = 0._WP

    G(2,1) = -V_2
    G(2,2) = V_2

    ! Finish

    return

  end function G_lagp_

  !****

  function H (this, s, x, omega)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: s
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: H(2,2)

    ! Evaluate the transformation matrix to convert variables to
    ! canonical form

    select case (this%set)
    case (DZIEM_SET)
       H = identity_matrix(2)
    case (JCD_SET)
       H = this%H_jcd_(s, x, omega)
    case (MIX_SET)
       H = identity_matrix(2)
    case (LAGP_SET)
       H = this%H_lagp_(s, x, omega)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function H

  !****

  function H_jcd_ (this, s, x, omega) result (H)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: s
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: H(2,2)

    real(WP) :: c_1
    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! to the canonical form

    ! Calculate coefficients

    c_1 = this%ml%c_1(s, x)

    omega_c = this%rt%omega_c(s, x, omega)

    ! Set up the matrix
      
    H(1,1) = 1._WP
    H(1,2) = 0._WP

    H(2,1) = 0._WP
    H(2,2) = c_1*omega_c**2

    ! Finish

    return

  end function H_jcd_

  !****

  function H_lagp_ (this, s, x, omega) result (H)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: s
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: H(2,2)

    real(WP) :: V_2

    ! Evaluate the transformation matrix to convert LAGP variables
    ! to the canonical form

    ! Calculate coefficients

    V_2 = this%ml%V_2(s, x)

    ! Set up the matrix

    H(1,1) = 1._WP
    H(1,2) = 0._WP

    H(2,1) = 1._WP
    H(2,2) = 1._WP/V_2

    ! Finish

    return

  end function H_lagp_

  !****

  function dH (this, s, x, omega)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: s
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: dH(2,2)

    ! Evaluate the derivative x dH/dx of the transformation matrix H

    select case (this%set)
    case (DZIEM_SET)
       dH = 0._WP
    case (JCD_SET)
       dH = this%dH_jcd_(s, x, omega)
    case (MIX_SET)
       dH = 0._WP
    case (LAGP_SET)
       dH = this%dH_lagp_(s, x, omega)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function dH

  !****

  function dH_jcd_ (this, s, x, omega) result (dH)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: s
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: dH(2,2)

    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: omega_c

    ! Evaluate the derivative x dH/dx of the JCD-variables
    ! transformation matrix H

    ! Calculate coefficients

    U = this%ml%U(s, x)
    c_1 = this%ml%c_1(s, x)

    omega_c = this%rt%omega_c(s, x, omega)

    ! Set up the matrix (nb: the derivative of omega_c is neglected;
    ! this is incorrect when rotation is non-zero)
      
    dH(1,1) = 0._WP
    dH(1,2) = 0._WP

    dH(2,1) = 0._WP
    dH(2,2) = c_1*(3._WP - U)*omega_c**2

    ! Finish

    return

  end function dH_jcd_

  !****

  function dH_lagp_ (this, s, x, omega) result (dH)

    class(rad_vars_t), intent(in) :: this
    integer, intent(in)           :: s
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: dH(2,2)

    real(WP) :: V_2
    real(WP) :: V
    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U

    ! Evaluate the derivative x dH/dx of the LAGP-variables
    ! transformation matrix H

    ! Calculate coefficients

    V_2 = this%ml%V_2(s, x)
    V = V_2*x**2
    V_g = V/this%ml%Gamma_1(s, x)
    As = this%ml%As(s, x)
    U = this%ml%U(s, x)

    ! Set up the matrix

    dH(1,1) = 0._WP
    dH(1,2) = 0._WP

    dH(2,1) = 0._WP
    dH(2,2) = -(-V_g - As + U + V - 3)/V_2

    ! Finish

    return

  end function dH_lagp_
  
end module gyre_rad_vars
