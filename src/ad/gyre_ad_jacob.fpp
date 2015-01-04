! Module   : gyre_ad_jacob
! Purpose  : jacobian evaluation (adiabatic)
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

module gyre_ad_jacob

  ! Uses

  use core_kinds

  use gyre_jacob
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

  ! Derived-type definitions

  type, extends (r_jacob_t) :: ad_jacob_t
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
     procedure, public :: T => T_
     procedure         :: T_jcd_
     procedure         :: T_mix_
  end type ad_jacob_t

  ! Interfaces

  interface ad_jacob_t
     module procedure ad_jacob_t_
  end interface ad_jacob_t

  ! Access specifiers

  private

  public :: ad_jacob_t

  ! Procedures

contains

  function ad_jacob_t_ (ml, rt, op) result (jc)

    class(model_t), pointer, intent(in) :: ml
    class(r_rot_t), intent(in)          :: rt
    type(osc_par_t), intent(in)         :: op
    type(ad_jacob_t)                    :: jc

    ! Construct the ad_jacob_t

    jc%ml => ml
    allocate(jc%rt, SOURCE=rt)

    select case (op%variables_set)
    case ('DZIEM')
       jc%vars = DZIEM_VARS
    case ('JCD')
       jc%vars = JCD_VARS
    case ('MIX')
       jc%vars = MIX_VARS
    case default
       $ABORT(Invalid variables_set)
    end select

    jc%n_e = 4

    ! Finish

    return

  end function ad_jacob_t_

!****

  function A_ (this, x, omega) result (A)

    class(ad_jacob_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: A(this%n_e,this%n_e)
    
    ! Evaluate the Jacobian matrix

    A = this%xA(x, omega)/x

    ! Finish

    return

  end function A_

!****

  function xA_ (this, x, omega) result (xA)

    class(ad_jacob_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: xA(this%n_e,this%n_e)
    
    ! Evaluate the log(x)-space Jacobian matrix (=x*A)

    select case (this%vars)
    case (DZIEM_VARS)
       xA = this%xA_dziem_(x, omega)
    case (JCD_VARS)
       xA = this%xA_jcd_(x, omega)
    case (MIX_VARS)
       xA = this%xA_mix_(x, omega)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function xA_

!****

  function xA_dziem_ (this, x, omega) result (xA)

    class(ad_jacob_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: xA(this%n_e,this%n_e)
    
    ! Evaluate the log(x)-space Jacobian matrix ([Dzi1971]
    ! formulation)

    associate(V_g => this%ml%V(x)/this%ml%Gamma_1(x), U => this%ml%U(x), &
              As => this%ml%As(x), c_1 => this%ml%c_1(x), &
              lambda => this%rt%lambda(x, omega), l_0 => this%rt%l_0(omega), omega_c => this%rt%omega_c(x, omega))

      xA(1,1) = V_g - 1._WP - l_0
      xA(1,2) = lambda/(c_1*omega_c**2) - V_g
      xA(1,3) = V_g
      xA(1,4) = 0._WP
      
      xA(2,1) = c_1*omega_c**2 - As
      xA(2,2) = As - U + 3._WP - l_0
      xA(2,3) = -As
      xA(2,4) = 0._WP
      
      xA(3,1) = 0._WP
      xA(3,2) = 0._WP
      xA(3,3) = 3._WP - U - l_0
      xA(3,4) = 1._WP
      
      xA(4,1) = U*As
      xA(4,2) = U*V_g
      xA(4,3) = lambda - U*V_g
      xA(4,4) = -U - l_0 + 2._WP

    end associate

    ! Finish

    return

  end function xA_dziem_

!****

  function xA_jcd_ (this, x, omega) result (xA)

    class(ad_jacob_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: xA(this%n_e,this%n_e)
    
    ! Evaluate the log(x)-space Jacobian matrix ([Chr2008]
    ! formulation)

    associate(V_g => this%ml%V(x)/this%ml%Gamma_1(x), U => this%ml%U(x), &
              As => this%ml%As(x), c_1 => this%ml%c_1(x), l => this%rt%mp%l, &
              lambda => this%rt%lambda(x, omega), l_0 => this%rt%l_0(omega), omega_c => this%rt%omega_c(x, omega))

      if (l /= 0) then

         xA(1,1) = V_g - 1._WP - l_0
         xA(1,2) = 1._WP - V_g*c_1*omega_c**2/lambda
         xA(1,3) = -V_g
         xA(1,4) = 0._WP
      
         xA(2,1) = lambda - As*lambda/(c_1*omega_c**2)
         xA(2,2) = As - l_0
         xA(2,3) = As*lambda/(c_1*omega_c**2)
         xA(2,4) = 0._WP
      
         xA(3,1) = 0._WP
         xA(3,2) = 0._WP
         xA(3,3) = 2._WP - l_0
         xA(3,4) = 1._WP
      
         xA(4,1) = -U*As
         xA(4,2) = -U*V_g*c_1*omega_c**2/lambda
         xA(4,3) = lambda + U*(As - 2._WP)
         xA(4,4) = 2._WP*(1._WP-U) - (l_0 - 1._WP)

      else

         xA(1,1) = V_g - 1._WP
         xA(1,2) = -V_g*c_1*omega_c**2
         xA(1,3) = -V_g
         xA(1,4) = 0._WP
      
         xA(2,1) = 1._WP - As/(c_1*omega_c**2)
         xA(2,2) = As
         xA(2,3) = As/(c_1*omega_c**2)
         xA(2,4) = 0._WP
      
         xA(3,1) = 0._WP
         xA(3,2) = 0._WP
         xA(3,3) = 2._WP
         xA(3,4) = 1._WP
      
         xA(4,1) = -U*As
         xA(4,2) = -U*V_g*c_1*omega_c**2
         xA(4,3) = U*(As - 2._WP)
         xA(4,4) = 2._WP*(1._WP-U) + 1._WP

      endif

    end associate

    ! Finish

    return

  end function xA_jcd_

!****

  function xA_mix_ (this, x, omega) result (xA)

    class(ad_jacob_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: xA(this%n_e,this%n_e)
    
    ! Evaluate the log(x)-space Jacobian matrix (mixed formulation)

    associate(V_g => this%ml%V(x)/this%ml%Gamma_1(x), U => this%ml%U(x), &
              As => this%ml%As(x), c_1 => this%ml%c_1(x), & 
              lambda => this%rt%lambda(x, omega), l_0 => this%rt%l_0(omega), omega_c => this%rt%omega_c(x, omega))

      xA(1,1) = V_g - 1._WP - l_0
      xA(1,2) = lambda/(c_1*omega_c**2) - V_g
      xA(1,3) = -V_g
      xA(1,4) = 0._WP
      
      xA(2,1) = c_1*omega_c**2 - As
      xA(2,2) = As - U + 3._WP - l_0
      xA(2,3) = As
      xA(2,4) = 0._WP
      
      xA(3,1) = 0._WP
      xA(3,2) = 0._WP
      xA(3,3) = 2._WP - l_0
      xA(3,4) = 1._WP
      
      xA(4,1) = -U*As
      xA(4,2) = -U*V_g
      xA(4,3) = lambda + U*(As - 2._WP)
      xA(4,4) = 2._WP*(1._WP-U) - (l_0 - 1._WP)

    end associate

    ! Finish

    return

  end function xA_mix_

!****

  function T_ (this, x, omega, to_canon) result (T)

    class(ad_jacob_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    logical, intent(in)           :: to_canon
    real(WP)                      :: T(this%n_e,this%n_e)

    ! Calculate the transformation matrix to convert variables to/from
    ! the canonical (DZEIM) formulation

    select case (this%vars)
    case (DZIEM_VARS)
       T = identity_matrix(this%n_e)
    case (JCD_VARS)
       T = this%T_jcd_(x, omega, to_canon)
    case (MIX_VARS)
       T = this%T_mix_(x, omega, to_canon)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function T_

!****

  function T_jcd_ (this, x, omega, to_canon) result (T)

    class(ad_jacob_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    logical, intent(in)           :: to_canon
    real(WP)                      :: T(this%n_e,this%n_e)

    ! Calculate the transformation matrix to convert JCD variables
    ! to/from the canonical (DZEIM) formulation

    associate(U => this%ml%U(x), c_1 => this%ml%c_1(x), l => this%rt%mp%l, &
              lambda => this%rt%lambda(x, omega), omega_c => this%rt%omega_c(x, omega))
      
      if (to_canon) then

         if (l /= 0._WP) then

            T(1,1) = 1._WP
            T(1,2) = 0._WP
            T(1,3) = 0._WP
            T(1,4) = 0._WP

            T(2,1) = 0._WP
            T(2,2) = c_1*omega_c**2/(lambda)
            T(2,3) = 0._WP
            T(2,4) = 0._WP

            T(3,1) = 0._WP
            T(3,2) = 0._WP
            T(3,3) = -1._WP
            T(3,4) = 0._WP

            T(4,1) = 0._WP
            T(4,2) = 0._WP
            T(4,3) = 1._WP - U
            T(4,4) = -1._WP

         else

            T(1,1) = 1._WP
            T(1,2) = 0._WP
            T(1,3) = 0._WP
            T(1,4) = 0._WP

            T(2,1) = 0._WP
            T(2,2) = c_1*omega_c**2
            T(2,3) = 0._WP
            T(2,4) = 0._WP

            T(3,1) = 0._WP
            T(3,2) = 0._WP
            T(3,3) = -1._WP
            T(3,4) = 0._WP

            T(4,1) = 0._WP
            T(4,2) = 0._WP
            T(4,3) = 1._WP - U
            T(4,4) = -1._WP

         endif

      else

         if (l /= 0) then

            T(1,1) = 1._WP
            T(1,2) = 0._WP
            T(1,3) = 0._WP
            T(1,4) = 0._WP

            T(2,1) = 0._WP
            T(2,2) = lambda/(c_1*omega_c**2)
            T(2,3) = 0._WP
            T(2,4) = 0._WP

            T(3,1) = 0._WP
            T(3,2) = 0._WP
            T(3,3) = -1._WP
            T(3,4) = 0._WP

            T(4,1) = 0._WP
            T(4,2) = 0._WP
            T(4,3) = -(1._WP - U)
            T(4,4) = -1._WP

         else

            T(1,1) = 1._WP
            T(1,2) = 0._WP
            T(1,3) = 0._WP
            T(1,4) = 0._WP

            T(2,1) = 0._WP
            T(2,2) = 1._WP/(c_1*omega_c**2)
            T(2,3) = 0._WP
            T(2,4) = 0._WP

            T(3,1) = 0._WP
            T(3,2) = 0._WP
            T(3,3) = -1._WP
            T(3,4) = 0._WP

            T(4,1) = 0._WP
            T(4,2) = 0._WP
            T(4,3) = -(1._WP - U)
            T(4,4) = -1._WP

         endif

      end if

    end associate

    ! Finish

    return

  end function T_jcd_

!****

  function T_mix_ (this, x, omega, to_canon) result (T)

    class(ad_jacob_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    logical, intent(in)           :: to_canon
    real(WP)                      :: T(this%n_e,this%n_e)

    ! Calculate the transformation matrix to convert MIX variables
    ! to/from the canonical (DZEIM) formulation

    associate(U => this%ml%U(x))

      if (to_canon) then

         T(1,1) = 1._WP
         T(1,2) = 0._WP
         T(1,3) = 0._WP
         T(1,4) = 0._WP

         T(2,1) = 0._WP
         T(2,2) = 1._WP
         T(2,3) = 0._WP
         T(2,4) = 0._WP

         T(3,1) = 0._WP
         T(3,2) = 0._WP
         T(3,3) = -1._WP
         T(3,4) = 0._WP

         T(4,1) = 0._WP
         T(4,2) = 0._WP
         T(4,3) = 1._WP - U
         T(4,4) = -1._WP

      else

         T(1,1) = 1._WP
         T(1,2) = 0._WP
         T(1,3) = 0._WP
         T(1,4) = 0._WP

         T(2,1) = 0._WP
         T(2,2) = 1._WP
         T(2,3) = 0._WP
         T(2,4) = 0._WP

         T(3,1) = 0._WP
         T(3,2) = 0._WP
         T(3,3) = -1._WP
         T(3,4) = 0._WP

         T(4,1) = 0._WP
         T(4,2) = 0._WP
         T(4,3) = -(1._WP - U)
         T(4,4) = -1._WP

      end if

    end associate

    ! Finish

    return

  end function T_mix_

end module gyre_ad_jacob
