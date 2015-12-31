! Module   : gyre_ad_vars
! Purpose  : adiabatic variables transformations
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

module gyre_ad_vars

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

  type :: ad_vars_t
     private
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     integer                     :: l
     integer                     :: set
   contains
     private
     procedure, public :: A => A_
     procedure         :: A_jcd_
     procedure         :: A_mix_
     procedure         :: A_lagp_
     procedure, public :: B => B_
     procedure         :: B_jcd_
     procedure         :: B_mix_
     procedure         :: B_lagp_
     procedure, public :: dB => dB_
     procedure         :: dB_jcd_
     procedure         :: dB_mix_
     procedure         :: dB_lagp_
  end type ad_vars_t

  ! Interfaces

  interface ad_vars_t
     module procedure ad_vars_t_
  end interface ad_vars_t

  ! Access specifiers

  private

  public :: ad_vars_t

  ! Procedures

contains

  function ad_vars_t_ (ml, md_p, os_p) result (vr)

    class(model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(ad_vars_t)                     :: vr

    ! Construct the ad_vars_t

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

    vr%l = md_p%l

    ! Finish

    return

  end function ad_vars_t_

  !****

  function A_ (this, s, x, omega) result (A)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: A(4,4)

    ! Evaluate the transformation matrix to convert variables from
    ! the canonical form

    select case (this%set)
    case (DZIEM_SET)
       A = identity_matrix(4)
    case (JCD_SET)
       A = this%A_jcd_(s, x, omega)
    case (MIX_SET)
       A = this%A_mix_(s, x, omega)
    case (LAGP_SET)
       A = this%A_lagp_(s, x, omega)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function A_

  !****

  function A_jcd_ (this, s, x, omega) result (A)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: A(4,4)

    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! from the canonical form

    ! Calculate coefficients

    U = this%ml%U(s, x)
    c_1 = this%ml%c_1(s, x)

    lambda = this%rt%lambda(s, x, omega)

    omega_c = this%rt%omega_c(s, x, omega)

    ! Set up the matrix
      
    if (this%l /= 0) then

       A(1,1) = 1._WP
       A(1,2) = 0._WP
       A(1,3) = 0._WP
       A(1,4) = 0._WP
          
       A(2,1) = 0._WP
       A(2,2) = lambda/(c_1*omega_c**2)
       A(2,3) = 0._WP
       A(2,4) = 0._WP
          
       A(3,1) = 0._WP
       A(3,2) = 0._WP
       A(3,3) = -1._WP
       A(3,4) = 0._WP

       A(4,1) = 0._WP
       A(4,2) = 0._WP
       A(4,3) = -(1._WP - U)
       A(4,4) = -1._WP

    else

       A(1,1) = 1._WP
       A(1,2) = 0._WP
       A(1,3) = 0._WP
       A(1,4) = 0._WP

       A(2,1) = 0._WP
       A(2,2) = 1._WP/(c_1*omega_c**2)
       A(2,3) = 0._WP
       A(2,4) = 0._WP

       A(3,1) = 0._WP
       A(3,2) = 0._WP
       A(3,3) = -1._WP
       A(3,4) = 0._WP

       A(4,1) = 0._WP
       A(4,2) = 0._WP
       A(4,3) = -(1._WP - U)
       A(4,4) = -1._WP

    endif

    ! Finish

    return

  end function A_jcd_

  !****

  function A_mix_ (this, s, x, omega) result (A)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: A(4,4)

    real(WP) :: U

    ! Evaluate the transformation matrix to convert MIX variables
    ! from the canonical form

    ! Calculate coefficients

    U = this%ml%U(s, x)

    ! Set up the matrix

    A(1,1) = 1._WP
    A(1,2) = 0._WP
    A(1,3) = 0._WP
    A(1,4) = 0._WP

    A(2,1) = 0._WP
    A(2,2) = 1._WP
    A(2,3) = 0._WP
    A(2,4) = 0._WP

    A(3,1) = 0._WP
    A(3,2) = 0._WP
    A(3,3) = -1._WP
    A(3,4) = 0._WP

    A(4,1) = 0._WP
    A(4,2) = 0._WP
    A(4,3) = -(1._WP - U)
    A(4,4) = -1._WP

    ! Finish

    return

  end function A_mix_

  !****

  function A_lagp_ (this, s, x, omega) result (A)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: A(4,4)

    real(WP) :: V_2

    ! Evaluate the transformation matrix to convert LAGP variables
    ! from the canonical form

    ! Calculate coefficients

    V_2 = this%ml%V_2(s, x)

    ! Set up the matrix

    A(1,1) = 1._WP
    A(1,2) = 0._WP
    A(1,3) = 0._WP
    A(1,4) = 0._WP

    A(2,1) = -V_2
    A(2,2) = V_2
    A(2,3) = -V_2
    A(2,4) = 0._WP

    A(3,1) = 0._WP
    A(3,2) = 0._WP
    A(3,3) = 1._WP
    A(3,4) = 0._WP

    A(4,1) = 0._WP
    A(4,2) = 0._WP
    A(4,3) = 0._WP
    A(4,4) = 1._WP

    ! Finish

    return

  end function A_lagp_

  !****

  function B_ (this, s, x, omega) result (B)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: B(4,4)

    ! Evaluate the transformation matrix to convert variables to
    ! canonical form

    select case (this%set)
    case (DZIEM_SET)
       B = identity_matrix(4)
    case (JCD_SET)
       B = this%B_jcd_(s, x, omega)
    case (MIX_SET)
       B = this%B_mix_(s, x, omega)
    case (LAGP_SET)
       B = this%B_lagp_(s, x, omega)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function B_

  !****

  function B_jcd_ (this, s, x, omega) result (B)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: B(4,4)

    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! to the canonical form

    ! Calculate coefficients

    U = this%ml%U(s, x)
    c_1 = this%ml%c_1(s, x)

    lambda = this%rt%lambda(s, x, omega)

    omega_c = this%rt%omega_c(s, x, omega)

    ! Set up the matrix
      
    if (this%l /= 0._WP) then

       B(1,1) = 1._WP
       B(1,2) = 0._WP
       B(1,3) = 0._WP
       B(1,4) = 0._WP
       
       B(2,1) = 0._WP
       B(2,2) = c_1*omega_c**2/lambda
       B(2,3) = 0._WP
       B(2,4) = 0._WP

       B(3,1) = 0._WP
       B(3,2) = 0._WP
       B(3,3) = -1._WP
       B(3,4) = 0._WP

       B(4,1) = 0._WP
       B(4,2) = 0._WP
       B(4,3) = 1._WP - U
       B(4,4) = -1._WP

    else

       B(1,1) = 1._WP
       B(1,2) = 0._WP
       B(1,3) = 0._WP
       B(1,4) = 0._WP

       B(2,1) = 0._WP
       B(2,2) = c_1*omega_c**2
       B(2,3) = 0._WP
       B(2,4) = 0._WP

       B(3,1) = 0._WP
       B(3,2) = 0._WP
       B(3,3) = -1._WP
       B(3,4) = 0._WP

       B(4,1) = 0._WP
       B(4,2) = 0._WP
       B(4,3) = 1._WP - U
       B(4,4) = -1._WP

    endif

    ! Finish

    return

  end function B_jcd_

  !****

  function B_mix_ (this, s, x, omega) result (B)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: B(4,4)

    real(WP) :: U

    ! Evaluate the transformation matrix to convert MIX variables
    ! to the canonical form

    ! Calculate coefficients

    U = this%ml%U(s, x)

    ! Set up the matrix

    B(1,1) = 1._WP
    B(1,2) = 0._WP
    B(1,3) = 0._WP
    B(1,4) = 0._WP

    B(2,1) = 0._WP
    B(2,2) = 1._WP
    B(2,3) = 0._WP
    B(2,4) = 0._WP

    B(3,1) = 0._WP
    B(3,2) = 0._WP
    B(3,3) = -1._WP
    B(3,4) = 0._WP

    B(4,1) = 0._WP
    B(4,2) = 0._WP
    B(4,3) = 1._WP - U
    B(4,4) = -1._WP

    ! Finish

    return

  end function B_mix_

  !****

  function B_lagp_ (this, s, x, omega) result (B)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: B(4,4)

    real(WP) :: V_2

    ! Evaluate the transformation matrix to convert LAGP variables
    ! to the canonical form

    ! Calculate coefficients

    V_2 = this%ml%V_2(s, x)

    ! Set up the matrix

    B(1,1) = 1._WP
    B(1,2) = 0._WP
    B(1,3) = 0._WP
    B(1,4) = 0._WP

    B(2,1) = 1._WP
    B(2,2) = 1._WP/V_2
    B(2,3) = 1._WP
    B(2,4) = 0._WP

    B(3,1) = 0._WP
    B(3,2) = 0._WP
    B(3,3) = 1._WP
    B(3,4) = 0._WP

    B(4,1) = 0._WP
    B(4,2) = 0._WP
    B(4,3) = 0._WP
    B(4,4) = 1._WP

    ! Finish

    return

  end function B_lagp_

  !****

  function dB_ (this, s, x, omega) result (dB)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: dB(4,4)

    ! Evaluate the derivative x dB/dx of the transformation matrix B

    select case (this%set)
    case (DZIEM_SET)
       dB = 0._WP
    case (JCD_SET)
       dB = this%dB_jcd_(s, x, omega)
    case (MIX_SET)
       dB = this%dB_mix_(s, x, omega)
    case (LAGP_SET)
       dB = this%dB_lagp_(s, x, omega)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function dB_

  !****

  function dB_jcd_ (this, s, x, omega) result (dB)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: dB(4,4)

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: omega_c

    ! Evaluate the derivative x dB/dx of the JCD-variables
    ! transformation matrix B

    ! Calculate coefficients

    V_g = this%ml%V_2(s, x)*x**2/this%ml%Gamma_1(s, x)
    As = this%ml%As(s, x) 
    U = this%ml%U(s, x)
    c_1 = this%ml%c_1(s, x)

    lambda = this%rt%lambda(s, x, omega)

    omega_c = this%rt%omega_c(s, x, omega)

    ! Set up the matrix (nb: the derivatives of omega_c and lambda are
    ! neglected; this is incorrect when rotation is non-zero)
      
    if (this%l /= 0._WP) then

       dB(1,1) = 0._WP
       dB(1,2) = 0._WP
       dB(1,3) = 0._WP
       dB(1,4) = 0._WP
       
       dB(2,1) = 0._WP
       dB(2,2) = c_1*(3._WP - U)*omega_c**2/lambda
       dB(2,3) = 0._WP
       dB(2,4) = 0._WP

       dB(3,1) = 0._WP
       dB(3,2) = 0._WP
       dB(3,3) = 0._WP
       dB(3,4) = 0._WP

       dB(4,1) = 0._WP
       dB(4,2) = 0._WP
       dB(4,3) = U*(V_g + As + U - 3._WP)
       dB(4,4) = 0._WP

    else

       dB(1,1) = 0._WP
       dB(1,2) = 0._WP
       dB(1,3) = 0._WP
       dB(1,4) = 0._WP

       dB(2,1) = 0._WP
       dB(2,2) = c_1*(3._WP - U)*omega_c**2
       dB(2,3) = 0._WP
       dB(2,4) = 0._WP

       dB(3,1) = 0._WP
       dB(3,2) = 0._WP
       dB(3,3) = 0._WP
       dB(3,4) = 0._WP

       dB(4,1) = 0._WP
       dB(4,2) = 0._WP
       dB(4,3) = U*(V_g + As + U - 3._WP)
       dB(4,4) = 0._WP

    endif

    ! Finish

    return

  end function dB_jcd_

  !****

  function dB_mix_ (this, s, x, omega) result (dB)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: dB(4,4)

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U

    ! Evaluate the derivative x dB/dx of the MIX-variables
    ! transformation matrix B

    ! Calculate coefficients

    V_g = this%ml%V_2(s, x)*x**2/this%ml%Gamma_1(s, x)
    As = this%ml%As(s, x)
    U = this%ml%U(s, x)

    ! Set up the matrix

    dB(1,1) = 0._WP
    dB(1,2) = 0._WP
    dB(1,3) = 0._WP
    dB(1,4) = 0._WP

    dB(2,1) = 0._WP
    dB(2,2) = 0._WP
    dB(2,3) = 0._WP
    dB(2,4) = 0._WP

    dB(3,1) = 0._WP
    dB(3,2) = 0._WP
    dB(3,3) = 0._WP
    dB(3,4) = 0._WP

    dB(4,1) = 0._WP
    dB(4,2) = 0._WP
    dB(4,3) = U*(V_g + As + U - 3._WP)
    dB(4,4) = 0._WP

    ! Finish

    return

  end function dB_mix_

!****

  function dB_lagp_ (this, s, x, omega) result (dB)

    class(ad_vars_t), intent(in) :: this
    integer, intent(in)          :: s
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: dB(4,4)

    real(WP) :: V_2
    real(WP) :: V
    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U

    ! Evaluate the derivative x dB/dx of the LAGP-variables
    ! transformation matrix B

    ! Calculate coefficients

    V_2 = this%ml%V_2(s, x)
    V = V_2*x**2
    V_g = V/this%ml%Gamma_1(s, x)
    As = this%ml%As(s, x)
    U = this%ml%U(s, x)

    ! Set up the matrix

    dB(1,1) = 0._WP
    dB(1,2) = 0._WP
    dB(1,3) = 0._WP
    dB(1,4) = 0._WP

    dB(2,1) = 0._WP
    dB(2,2) = -(-V_g - As + U + V - 3)/V_2
    dB(2,3) = 0._WP
    dB(2,4) = 0._WP

    dB(3,1) = 0._WP
    dB(3,2) = 0._WP
    dB(3,3) = 0._WP
    dB(3,4) = 0._WP

    dB(4,1) = 0._WP
    dB(4,2) = 0._WP
    dB(4,3) = 0._WP
    dB(4,4) = 0._WP

    ! Finish

    return

  end function dB_lagp_

end module gyre_ad_vars
