! Module   : gyre_nad_bound
! Purpose  : nonadiabatic boundary conditions
!
! Copyright 2013 Rich Townsend
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

module gyre_nad_bound

  ! Uses

  use core_kinds

  use gyre_base_coeffs
  use gyre_therm_coeffs
  use gyre_oscpar
  use gyre_ad_bound

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: nad_bound_t
     private
     class(base_coeffs_t), pointer  :: bc => null()
     class(therm_coeffs_t), pointer :: tc => null()
     type(oscpar_t), pointer        :: op => null()
     integer, public                :: n_e
     integer, public                :: n_i
     integer, public                :: n_o
   contains 
     private
     procedure, public :: init
     procedure, public :: inner_bound
     procedure, public :: outer_bound
     procedure, public :: outer_bound_zero
     procedure, public :: outer_bound_dziem
     procedure, public :: outer_bound_unno
     procedure, public :: outer_bound_jcd
  end type nad_bound_t

  ! Access specifiers

  private

  public :: nad_bound_t

  ! Procedures

contains

  subroutine init (this, bc, tc, op)

    class(nad_bound_t), intent(out)           :: this
    class(base_coeffs_t), intent(in), target  :: bc
    class(therm_coeffs_t), intent(in), target :: tc
    type(oscpar_t), intent(in), target        :: op

    ! Initialize the nad_bound

    this%bc => bc
    this%tc => tc
    this%op => op

    this%n_i = 3
    this%n_o = 3
    this%n_e = this%n_i + this%n_o

    ! Finish

    return
    
  end subroutine init

!****

  function inner_bound (this, x_i, omega) result (B_i)

    class(nad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: x_i
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_i(this%n_i,this%n_e)

    $ASSERT(x_i == 0._WP,Boundary condition invalid for x_i /= 0)

    ! Set the inner boundary conditions to enforce non-diverging modes

    associate(c_1 => this%bc%c_1(x_i), l => this%op%l)

      B_i(1,1) = c_1*omega**2
      B_i(1,2) = -l
      B_i(1,3) = 0._WP
      B_i(1,4) = 0._WP
      B_i(1,5) = 0._WP
      B_i(1,6) = 0._WP

      B_i(2,1) = 0._WP
      B_i(2,2) = 0._WP
      B_i(2,3) = l
      B_i(2,4) = -1._WP
      B_i(2,5) = 0._WP
      B_i(2,6) = 0._WP

      B_i(3,1) = 0._WP
      B_i(3,2) = 0._WP
      B_i(3,3) = 0._WP
      B_i(3,4) = 0._WP
      B_i(3,5) = 1._WP
      B_i(3,6) = 0._WP

    end associate

    ! Finish

    return

  end function inner_bound

!****

  function outer_bound (this, x_o, omega) result (B_o)

    class(nad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: x_o
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    ! Set the outer boundary conditions

    select case (this%op%outer_bound_type)
    case ('ZERO')
       B_o = this%outer_bound_zero(x_o, omega)
    case ('DZIEM')
       B_o = this%outer_bound_dziem(x_o, omega)
    case ('UNNO')
       B_o = this%outer_bound_unno(x_o, omega)
    case ('JCD')
       B_o = this%outer_bound_jcd(x_o, omega)
    case default
       $ABORT(Invalid outer_bound_type)
    end select

    ! Finish

    return

  end function outer_bound

!****

  function outer_bound_zero (this, x_o, omega) result (B_o)

    class(nad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: x_o
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    ! Set the outer boundary conditions, assuming delta p -> 0. The U
    ! term in the gravitational bc is required for cases where the
    ! surface density remains finite (see Cox 1980, eqn. 17.71)

    associate(V => this%bc%V(x_o), U => this%bc%U(x_o), nabla_ad => this%bc%nabla_ad(x_o), &
              l => this%op%l)

      B_o(1,1) = 1._WP
      B_o(1,2) = -1._WP
      B_o(1,3) = 1._WP
      B_o(1,4) = 0._WP
      B_o(1,5) = 0._WP
      B_o(1,6) = 0._WP
      
      B_o(2,1) = U
      B_o(2,2) = 0._WP
      B_o(2,3) = l + 1._WP
      B_o(2,4) = 1._WP
      B_o(2,5) = 0._WP
      B_o(2,6) = 0._WP

      B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
      B_o(3,2) = 4._WP*nabla_ad*V
      B_o(3,3) = -4._WP*nabla_ad*V
      B_o(3,4) = 0._WP
      B_o(3,5) = 4._WP
      B_o(3,6) = -1._WP

    end associate

    ! Finish

    return

  end function outer_bound_zero

!****

  function outer_bound_dziem (this, x_o, omega) result (B_o)

    class(nad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: x_o
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    ! Set the outer boundary conditions, assuming Dziembowski's (1971)
    ! condition: d(delta p)/dr -> 0 for an isothermal atmosphere.

    associate(V => this%bc%V(x_o), nabla_ad => this%bc%nabla_ad(x_o), &
              l => this%op%l)

      B_o(1,1) = 1 + (l*(l+1)/omega**2 - 4 - omega**2)/V
      B_o(1,2) = -1._WP
      B_o(1,3) = 1 + (l*(l+1)/omega**2 - l - 1)/V
      B_o(1,4) = 0._WP
      B_o(1,5) = 0._WP
      B_o(1,6) = 0._WP
     
      B_o(2,1) = 0._WP
      B_o(2,2) = 0._WP
      B_o(2,3) = l + 1._WP
      B_o(2,4) = 1._WP
      B_o(2,5) = 0._WP
      B_o(2,6) = 0._WP

      B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
      B_o(3,2) = 4._WP*nabla_ad*V
      B_o(3,3) = -4._WP*nabla_ad*V
      B_o(3,4) = 0._WP
      B_o(3,5) = 4._WP
      B_o(3,6) = -1._WP

    end associate

    ! Finish

    return

  end function outer_bound_dziem

!****

  function outer_bound_unno (this, x_o, omega) result (B_o)

    class(nad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: x_o
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    real(WP)    :: V_g
    real(WP)    :: As
    real(WP)    :: c_1
    complex(WP) :: lambda
    complex(WP) :: b_11
    complex(WP) :: b_12
    complex(WP) :: b_13
    complex(WP) :: b_21
    complex(WP) :: b_22
    complex(WP) :: b_23
    complex(WP) :: alpha_1
    complex(WP) :: alpha_2

    ! Set the outer boundary conditions, assuming Unno et al.'s (1989,
    ! S18.1) formulation.

    call eval_outer_coeffs_unno(this%bc, x_o, V_g, As, c_1)

    associate(V => this%bc%V(x_o), nabla_ad => this%bc%nabla_ad(x_o), l => this%op%l)

      lambda = outer_wavenumber(V_g, As, c_1, omega, l)
      
      b_11 = V_g - 3._WP
      b_12 = l*(l+1)/(c_1*omega**2) - V_g
      b_13 = V_g

      b_21 = c_1*omega**2 - As
      b_22 = 1._WP + As
      b_23 = -As
    
      alpha_1 = (b_12*b_23 - b_13*(b_22+l))/((b_11+l)*(b_22+l) - b_12*b_21)
      alpha_2 = (b_21*b_13 - b_23*(b_11+l))/((b_11+l)*(b_22+l) - b_12*b_21)

      B_o(1,1) = (lambda - b_11)/b_12
      B_o(1,2) = -1._WP
      B_o(1,3) = -(alpha_1*(lambda - b_11)/b_12 - alpha_2)
      B_o(1,4) = 0._WP
      B_o(1,5) = 0._WP
      B_o(1,6) = 0._WP

      B_o(2,1) = 0._WP
      B_o(2,2) = 0._WP
      B_o(2,3) = l + 1._WP
      B_o(2,4) = 1._WP
      B_o(2,5) = 0._WP
      B_o(2,6) = 0._WP

      B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
      B_o(3,2) = 4._WP*nabla_ad*V
      B_o(3,3) = -4._WP*nabla_ad*V
      B_o(3,4) = 0._WP
      B_o(3,5) = 4._WP
      B_o(3,6) = -1._WP

    end associate

    ! Finish

    return

  end function outer_bound_unno

!****

  function outer_bound_jcd (this, x_o, omega) result (B_o)

    class(nad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: x_o
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    real(WP)    :: V_g
    real(WP)    :: As
    real(WP)    :: c_1
    complex(WP) :: lambda
    complex(WP) :: b_11
    complex(WP) :: b_12

    ! Set the outer boundary conditions, assuming
    ! Christensen-Dalsgaard's formulation (see ADIPLS documentation)

    call eval_outer_coeffs_jcd(this%bc, x_o, V_g, As, c_1)

    associate(V => this%bc%V(x_o), nabla_ad => this%bc%nabla_ad(x_o), l => this%op%l)

      lambda = outer_wavenumber(V_g, As, c_1, omega, l)

      b_11 = V_g - 3._WP
      b_12 = l*(l+1)/(c_1*omega**2) - V_g

      if(l /= 0) then
         B_o(1,1) = (lambda - b_11)/b_12
         B_o(1,2) = -1._WP
         B_o(1,3) = 1._WP + (l*(l+1)/(c_1*omega**2) - l - 1._WP)/(V_g + As)
         B_o(1,4) = 0._WP
         B_o(1,5) = 0._WP
         B_o(1,6) = 0._WP
      else
         B_o(1,1) = (lambda - b_11)/b_12
         B_o(1,2) = -1._WP
         B_o(1,3) = 1._WP
         B_o(1,4) = 0._WP
         B_o(1,5) = 0._WP
         B_o(1,6) = 0._WP
      endif

      B_o(2,1) = 0._WP
      B_o(2,2) = 0._WP
      B_o(2,3) = l + 1._WP
      B_o(2,4) = 1._WP
      B_o(2,5) = 0._WP
      B_o(2,6) = 0._WP

      B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
      B_o(3,2) = 4._WP*nabla_ad*V
      B_o(3,3) = -4._WP*nabla_ad*V
      B_o(3,4) = 0._WP
      B_o(3,5) = 4._WP
      B_o(3,6) = -1._WP

    end associate

    ! Finish

    return

  end function outer_bound_jcd

end module gyre_nad_bound
