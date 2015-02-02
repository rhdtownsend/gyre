! Module   : gyre_rad_eqns
! Purpose  : differential equations evaluation (adiabatic radial)
!
! Copyright 2013-2015 Rich Townsend
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

module gyre_rad_eqns

  ! Uses

  use core_kinds

  use gyre_eqns
  use gyre_linalg
  use gyre_model
  use gyre_osc_par
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: DZIEM_VARS = 1
  integer, parameter :: JCD_VARS = 2
  integer, parameter :: MIX_VARS = 3
  integer, parameter :: LAGP_VARS = 4

  ! Derived-type definitions

  type, extends (r_eqns_t) :: rad_eqns_t
     private
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     integer                     :: vars
   contains
     private
     procedure, public :: A => A_
     procedure, public :: xA => xA_
     procedure         :: xA_dziem_
     procedure         :: xA_jcd_
     procedure         :: xA_mix_
     procedure         :: xA_lagp_
     procedure, public :: T => T_
     procedure         :: T_jcd_
     procedure         :: T_lagp_
  end type rad_eqns_t

  ! Interfaces

  interface rad_eqns_t
     module procedure rad_eqns_t_
  end interface rad_eqns_t

  ! Access specifiers

  private

  public :: rad_eqns_t

  ! Procedures

contains

  function rad_eqns_t_ (ml, rt, op) result (eq)

    class(model_t), pointer, intent(in) :: ml
    class(r_rot_t), intent(in)          :: rt
    type(osc_par_t), intent(in)         :: op
    type(rad_eqns_t)                    :: eq

    ! Construct the rad_eqns_t

    eq%ml => ml
    allocate(eq%rt, SOURCE=rt)

    select case (op%variables_set)
    case ('DZIEM')
       eq%vars = DZIEM_VARS
    case ('JCD')
       eq%vars = JCD_VARS
    case ('MIX')
       eq%vars = MIX_VARS
    case ('LAGP')
       eq%vars = LAGP_VARS
    case default
       $ABORT(Invalid variables_set)
    end select

    eq%n_e = 2

    ! Finish

    return

  end function rad_eqns_t_

!****

  function A_ (this, x, omega) result (A)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(x, omega)/x

    ! Finish

    return

  end function A_

!****

  function xA_ (this, x, omega) result (xA)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: xA(this%n_e,this%n_e)
    
    ! Evaluate the log(x)-space RHS matrix (=x*A)

    select case (this%vars)
    case (DZIEM_VARS)
       xA = this%xA_dziem_(x, omega)
    case (JCD_VARS)
       xA = this%xA_jcd_(x, omega)
    case (MIX_VARS)
       xA = this%xA_mix_(x, omega)
    case (LAGP_VARS)
       xA = this%xA_lagp_(x, omega)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function xA_

!****

  function xA_dziem_ (this, x, omega) result (xA)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: xA(this%n_e,this%n_e)

    real(WP) :: V_g
    real(WP) :: U
    real(WP) :: As
    real(WP) :: c_1
    real(WP) :: omega_c
    
    ! Evaluate the log(x)-space RHS matrix ([Dzi1971] formulation)

    ! Calculate coefficients

    V_g = this%ml%V_2(x)*x**2/this%ml%Gamma_1(x)
    U = this%ml%U(x)
    As = this%ml%As(x)
    c_1 = this%ml%c_1(x)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix

    xA(1,1) = V_g - 1._WP
    xA(1,2) = -V_g
      
    xA(2,1) = c_1*omega_c**2 + U - As
    xA(2,2) = As - U + 3._WP

    ! Finish

    return

  end function xA_dziem_

!****

  function xA_jcd_ (this, x, omega) result (xA)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: xA(this%n_e,this%n_e)
    
    real(WP) :: V_g
    real(WP) :: U
    real(WP) :: As
    real(WP) :: c_1
    real(WP) :: omega_c
    
    ! Evaluate the log(x)-space RHS matrix ([ChrDal2008] formulation)

    ! Calculate coefficients

    V_g = this%ml%V_2(x)*x**2/this%ml%Gamma_1(x)
    U = this%ml%U(x)
    As = this%ml%As(x)
    c_1 = this%ml%c_1(x)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix

    xA(1,1) = V_g - 1._WP
    xA(1,2) = -V_g*c_1*omega_c**2
    
    xA(2,1) = 1._WP - (As - U)/(c_1*omega_c**2)
    xA(2,2) = As
    
    ! Finish

    return

  end function xA_jcd_

!****

  function xA_mix_ (this, x, omega) result (xA)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: xA(this%n_e,this%n_e)
    
    real(WP) :: V_g
    real(WP) :: U
    real(WP) :: As
    real(WP) :: c_1
    real(WP) :: omega_c
    
    ! Evaluate the log(x)-space RHS matrix (mixed formulation)

    ! Calculate coefficients

    V_g = this%ml%V_2(x)*x**2/this%ml%Gamma_1(x)
    U = this%ml%U(x)
    As = this%ml%As(x)
    c_1 = this%ml%c_1(x)

    omega_c = this%rt%omega_c(x, omega)

    ! Set upt the matrix

    xA(1,1) = V_g - 1._WP
    xA(1,2) = -V_g
      
    xA(2,1) = c_1*omega_c**2 + U - As
    xA(2,2) = As - U + 3._WP

    ! Finish

    return

  end function xA_mix_

!****

  function xA_lagp_ (this, x, omega) result (xA)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: xA(this%n_e,this%n_e)
    
    real(WP) :: V_2
    real(WP) :: Gamma_1
    real(WP) :: c_1
    real(WP) :: omega_c

    ! Evaluate the log(x)-space RHS matrix (Lagrangian pressure
    ! perturbation formulation)

    ! Calculate coefficients

    V_2 = this%ml%V_2(x)
    Gamma_1 = this%ml%Gamma_1(x)
    c_1 = this%ml%c_1(x)

    omega_c = this%rt%omega_c(x, omega)

    xA(1,1) = -1._WP
    xA(1,2) = -x**2/Gamma_1
      
    xA(2,1) = V_2*(4._WP + c_1*omega_c**2)
    xA(2,2) = V_2*x**2

    ! Finish

    return

  end function xA_lagp_

!****

  function T_ (this, x, omega, to_canon) result (T)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    logical, intent(in)           :: to_canon
    real(WP)                      :: T(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert variables to/from
    ! the canonical (DZEIM) formulation

    select case (this%vars)
    case (DZIEM_VARS)
       T = identity_matrix(this%n_e)
    case (JCD_VARS)
       T = this%T_jcd_(x, omega, to_canon)
    case (MIX_VARS)
       T = identity_matrix(this%n_e)
    case (LAGP_VARS)
       T = this%T_lagp_(x, omega, to_canon)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function T_

!****

  function T_jcd_ (this, x, omega, to_canon) result (T)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    logical, intent(in)           :: to_canon
    real(WP)                      :: T(this%n_e,this%n_e)

    real(WP) :: c_1
    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! to/from the canonical (DZEIM) formulation

    ! Calculate coefficients

    c_1 = this%ml%c_1(x)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix

    if (to_canon) then

       T(1,1) = 1._WP
       T(1,2) = 0._WP
       
       T(2,1) = 0._WP
       T(2,2) = c_1*omega_c**2

    else

       T(1,1) = 1._WP
       T(1,2) = 0._WP

       T(2,1) = 0._WP
       T(2,2) = 1._WP/(c_1*omega_c**2)

    end if

    ! Finish

    return

  end function T_jcd_

!****

  function T_lagp_ (this, x, omega, to_canon) result (T)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    logical, intent(in)           :: to_canon
    real(WP)                      :: T(this%n_e,this%n_e)

    real(WP) :: V_2

    ! Evaluate the transformation matrix to convert LAGP variables
    ! to/from the canonical (DZEIM) formulation

    ! Calculate coefficients

    V_2 = this%ml%V_2(x)

    ! Set up the matrix

    if (to_canon) then

       T(1,1) = 1._WP
       T(1,2) = 0._WP

       T(2,1) = 1._WP
       T(2,2) = 1._WP/V_2

    else

       T(1,1) = 1._WP
       T(1,2) = 0._WP

       T(2,1) = -V_2
       T(2,2) = V_2

    end if

    ! Finish

    return

  end function T_lagp_

end module gyre_rad_eqns
