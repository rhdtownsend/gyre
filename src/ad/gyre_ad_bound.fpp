! Incfile  : gyre_ad_bound
! Purpose  : boundary conditions (adiabatic)
!
! Copyright 2013-2014 Rich Townsend
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

module gyre_ad_bound

  ! Uses

  use core_kinds

  use gyre_atmos
  use gyre_bound
  use gyre_jacob
  use gyre_model
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: REGULAR_TYPE_I = 1
  integer, parameter :: ZERO_TYPE_I = 2

  integer, parameter :: ZERO_TYPE_O = 3
  integer, parameter :: DZIEM_TYPE_O = 4
  integer, parameter :: UNNO_TYPE_O = 5
  integer, parameter :: JCD_TYPE_O = 6

  ! Derived-type definitions

  type, extends (r_bound_t) :: ad_bound_t
     private
     class(model_t), pointer       :: ml => null()
     class(r_rot_t), allocatable   :: rt
     class(r_jacob_t), allocatable :: jc
     real(WP)                      :: x_i
     real(WP)                      :: x_o
     integer                       :: type_i
     integer                       :: type_o
   contains 
     private
     procedure, public :: B_i => B_i_
     procedure         :: B_i_regular_
     procedure         :: B_i_zero_
     procedure, public :: B_o => B_o_
     procedure         :: B_o_zero_
     procedure         :: B_o_dziem_
     procedure         :: B_o_unno_
     procedure         :: B_o_jcd_
  end type ad_bound_t

  ! Interfaces

  interface ad_bound_t
     module procedure ad_bound_t_
  end interface ad_bound_t

  ! Access specifiers

  private

  public :: ad_bound_t

  ! Procedures

contains

  function ad_bound_t_ (ml, rt, jc, x_i, x_o, type_i, type_o) result (bd)

    class(model_t), pointer, intent(in) :: ml
    class(r_rot_t), intent(in)          :: rt
    class(r_jacob_t), intent(in)        :: jc
    real(WP)                            :: x_i
    real(WP)                            :: x_o
    character(*), intent(in)            :: type_i
    character(*), intent(in)            :: type_o
    type(ad_bound_t)                    :: bd

    ! Construct the ad_bound_t

    bd%ml => ml
    allocate(bd%rt, SOURCE=rt)
    allocate(bd%jc, SOURCE=jc)

    bd%x_i = x_i
    bd%x_o = x_o

    select case (type_i)
    case ('REGULAR')
       bd%type_i = REGULAR_TYPE_I
    case ('ZERO')
       bd%type_i = ZERO_TYPE_I
    case default
       $ABORT(Invalid type_i)
    end select

    select case (type_o)
    case ('ZERO')
       bd%type_o = ZERO_TYPE_O
    case ('DZIEM')
       bd%type_o = DZIEM_TYPE_O
    case ('UNNO')
       bd%type_o = UNNO_TYPE_O
    case ('JCD')
       bd%type_o = JCD_TYPE_O
    case default
       $ABORT(Invalid type_o)
    end select

    bd%n_i = 2
    bd%n_o = 2
    bd%n_e = bd%n_i + bd%n_o

    ! Finish

    return
    
  end function ad_bound_t_

!****

  function B_i_ (this, omega) result (B_i)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: omega
    real(WP)                      :: B_i(this%n_i,this%n_e)

    ! Evaluate the inner boundary conditions

    select case (this%type_i)
    case (REGULAR_TYPE_I)
       B_i = this%B_i_regular_(omega)
    case (ZERO_TYPE_I)
       B_i = this%B_i_zero_(omega)
    case default
       $ABORT(Invalid type_i)
    end select

    ! Transform to the variables used in the jacobian

    B_i = MATMUL(B_i, this%jc%T(this%x_i, omega, .TRUE.))

    ! Finish

    return

  end function B_i_

!****

  function B_i_regular_ (this, omega) result (B_i)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: omega
    real(WP)                      :: B_i(this%n_i,this%n_e)

    $ASSERT(this%x_i == 0._WP,Boundary condition invalid for x_i /= 0)

    ! Evaluate the inner boundary conditions (regular-enforcing)

    associate(c_1 => this%ml%c_1(this%x_i), &
              l_e => this%rt%l_e(this%x_i, omega), omega_c => this%rt%omega_c(this%x_i, omega))
                 
      B_i(1,1) = c_1*omega_c**2
      B_i(1,2) = -l_e
      B_i(1,3) = 0._WP
      B_i(1,4) = 0._WP
        
      B_i(2,1) = 0._WP
      B_i(2,2) = 0._WP
      B_i(2,3) = l_e
      B_i(2,4) = -1._WP
      
    end associate

    ! Finish

    return

  end function B_i_regular_

!****

  function B_i_zero_ (this, omega) result (B_i)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: omega
    real(WP)                      :: B_i(this%n_i,this%n_e)

    $ASSERT(this%x_i /= 0._WP,Boundary condition invalid for x_i == 0)

    ! Evaluate the inner boundary conditions (zero displacement/gravity)

    B_i(1,1) = 1._WP
    B_i(1,2) = 0._WP
    B_i(1,3) = 0._WP
    B_i(1,4) = 0._WP
        
    B_i(2,1) = 0._WP
    B_i(2,2) = 0._WP
    B_i(2,3) = 0._WP
    B_i(2,4) = 1._WP
      
    ! Finish

    return

  end function B_i_zero_

!****

  function B_o_ (this, omega) result (B_o)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: omega
    real(WP)                      :: B_o(this%n_o,this%n_e)

    ! Evaluate the outer boundary conditions

    select case (this%type_o)
    case (ZERO_TYPE_O)
       B_o = this%B_o_zero_(omega)
    case (DZIEM_TYPE_O)
       B_o = this%B_o_dziem_(omega)
    case (UNNO_TYPE_O)
       B_o = this%B_o_unno_(omega)
    case (JCD_TYPE_O)
       B_o = this%B_o_jcd_(omega)
    case default
       $ABORT(Invalid type_o)
    end select

    ! Transform to the variables used in the jacobian

    B_o = MATMUL(B_o, this%jc%T(this%x_o, omega, .TRUE.))
    
    ! Finish

    return

  end function B_o_

!****

  function B_o_zero_ (this, omega) result (B_o)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: omega
    real(WP)                      :: B_o(this%n_o,this%n_e)

    ! Evaluate the outer boundary conditions (zero-pressure)

    associate(U => this%ml%U(this%x_o), &
              l_e => this%rt%l_e(this%x_o, omega))

      B_o(1,1) = 1._WP
      B_o(1,2) = -1._WP
      B_o(1,3) = 1._WP
      B_o(1,4) = 0._WP
      
      B_o(2,1) = U
      B_o(2,2) = 0._WP
      B_o(2,3) = l_e + 1._WP
      B_o(2,4) = 1._WP

    end associate

    ! Finish

    return

  end function B_o_zero_

!****

  function B_o_dziem_ (this, omega) result (B_o)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: omega
    real(WP)                      :: B_o(this%n_o,this%n_e)

    ! Evaluate the outer boundary conditions ([Dzi1971] formulation)

    associate(V => this%ml%V(this%x_o), c_1 => this%ml%c_1(this%x_o), &
              l_e => this%rt%l_e(this%x_o, omega), omega_c => this%rt%omega_c(this%x_o, omega))
        
      B_o(1,1) = 1 + (l_e*(l_e+1._WP)/(c_1*omega_c**2) - 4._WP - c_1*omega_c**2)/V
      B_o(1,2) = -1._WP
      B_o(1,3) = 1 + (l_e*(l_e+1._WP)/(c_1*omega_c**2) - l_e - 1._WP)/V
      B_o(1,4) = 0._WP
      
      B_o(2,1) = 0._WP
      B_o(2,2) = 0._WP
      B_o(2,3) = l_e + 1._WP
      B_o(2,4) = 1._WP

    end associate

    ! Finish

    return

  end function B_o_dziem_

!****

  function B_o_unno_ (this, omega) result (B_o)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: omega
    real(WP)                      :: B_o(this%n_o,this%n_e)

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: b_11
    real(WP) :: b_12
    real(WP) :: b_13
    real(WP) :: b_21
    real(WP) :: b_22
    real(WP) :: b_23
    real(WP) :: alpha_1
    real(WP) :: alpha_2

    ! Evaluate the outer boundary conditions ([Unn1989] formulation)

    call eval_atmos_coeffs_unno(this%ml, this%x_o, V_g, As, c_1)

    associate(l_e => this%rt%l_e(this%x_o, omega), omega_c => this%rt%omega_c(this%x_o, omega))

      lambda = atmos_wavenumber(V_g, As, c_1, omega_c, l_e)
      
      b_11 = V_g - 3._WP
      b_12 = l_e*(l_e+1._WP)/(c_1*omega_c**2) - V_g
      b_13 = V_g

      b_21 = c_1*omega_c**2 - As
      b_22 = 1._WP + As
      b_23 = -As
    
      alpha_1 = (b_12*b_23 - b_13*(b_22+l_e))/((b_11+l_e)*(b_22+l_e) - b_12*b_21)
      alpha_2 = (b_21*b_13 - b_23*(b_11+l_e))/((b_11+l_e)*(b_22+l_e) - b_12*b_21)

      B_o(1,1) = lambda - b_11
      B_o(1,2) = -b_12
      B_o(1,3) = -(alpha_1*(lambda - b_11) - alpha_2*b_12)
      B_o(1,4) = 0._WP

      B_o(2,1) = 0._WP
      B_o(2,2) = 0._WP
      B_o(2,3) = l_e + 1._WP
      B_o(2,4) = 1._WP

    end associate

    ! Finish

    return

  end function B_o_unno_

!****

  function B_o_jcd_ (this, omega) result (B_o)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: omega
    real(WP)                      :: B_o(this%n_o,this%n_e)

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: b_11
    real(WP) :: b_12

    ! Evaluate the outer boundary conditions ([Chr2008] formulation)

    call eval_atmos_coeffs_jcd(this%ml, this%x_o, V_g, As, c_1)

    associate(l_e => this%rt%l_e(this%x_o, omega), omega_c => this%rt%omega_c(this%x_o, omega))

      lambda = atmos_wavenumber(V_g, As, c_1, omega_c, l_e)

      b_11 = V_g - 3._WP
      b_12 = l_e*(l_e+1._WP)/(c_1*omega_c**2) - V_g

      B_o(1,1) = lambda - b_11
      B_o(1,2) = -b_12
      B_o(1,3) = b_12 + (l_e*(l_e+1._WP)/(c_1*omega_c**2) - l_e - 1._WP)*b_12/(V_g + As)
      B_o(1,4) = 0._WP

      B_o(2,1) = 0._WP
      B_o(2,2) = 0._WP
      B_o(2,3) = l_e + 1._WP
      B_o(2,4) = 1._WP

    end associate

    ! Finish

    return

  end function B_o_jcd_

end module gyre_ad_bound
