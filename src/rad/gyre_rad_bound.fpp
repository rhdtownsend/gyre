! Module   : gyre_rad_bound
! Purpose  : radial adiabatic boundary conditions
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

module gyre_rad_bound

  ! Uses

  use core_kinds

  use gyre_base_coeffs
  use gyre_oscpar
  use gyre_ad_bound

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: rad_bound_t
     private
     class(base_coeffs_t), pointer :: bc => null()
     type(oscpar_t), pointer       :: op => null()
     integer, public               :: n_e
     integer, public               :: n_i
     integer, public               :: n_o
   contains 
     private
     procedure, public :: init
     procedure, public :: inner_bound
     procedure, public :: outer_bound
     procedure, public :: outer_bound_zero
     procedure, public :: outer_bound_dziem
     procedure, public :: outer_bound_unno
     procedure, public :: outer_bound_jcd
  end type rad_bound_t
  
  ! Access specifiers

  private

  public :: rad_bound_t

  ! Procedures

contains

  subroutine init (this, bc, op)

    class(rad_bound_t), intent(out)          :: this
    class(base_coeffs_t), intent(in), target :: bc
    type(oscpar_t), intent(in), target       :: op

    ! Initialize the ad_bound

    this%bc => bc
    this%op => op

    this%n_i = 1
    this%n_o = 1
    this%n_e = this%n_i + this%n_o

    ! Finish

    return
    
  end subroutine init

!****

  function inner_bound (this, x_i, omega) result (B_i)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: x_i
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_i(this%n_i,this%n_e)

    $ASSERT(x_i == 0._WP,Boundary condition invalid for x_i /= 0)

    ! Set the inner boundary conditions to enforce non-diverging modes

    B_i(1,1) = 1._WP
    B_i(1,2) = 0._WP

    ! Finish

    return

  end function inner_bound

!****

  function outer_bound (this, x_o, omega) result (B_o)

    class(rad_bound_t), intent(in) :: this
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

    class(rad_bound_t), intent(in) :: this
    real(WP)                       :: x_o
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    ! Set the outer boundary conditions, assuming delta p -> 0

    B_o(1,1) = 1._WP
    B_o(1,2) = -1._WP

    ! Finish

    return

  end function outer_bound_zero

!****

  function outer_bound_dziem (this, x_o, omega) result (B_o)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: x_o
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    ! Set the outer boundary conditions, assuming Dziembowski's (1971)
    ! condition: d(delta p)/dr -> 0 for an isothermal atmosphere

    associate(V => this%bc%V(x_o), c_1 => this%bc%V(x_o), &
              omega_c => omega)
        
      B_o(1,1) = 1 - (4._WP + c_1*omega_c**2)/V
      B_o(1,2) = -1._WP

    end associate

    ! Finish

    return

  end function outer_bound_dziem

!****

  function outer_bound_unno (this, x_o, omega) result (B_o)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: x_o
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    real(WP)    :: V_g
    real(WP)    :: As
    real(WP)    :: c_1
    complex(WP) :: lambda
    complex(WP) :: b_11
    complex(WP) :: b_12

    ! Set the outer boundary conditions, assuming Unno et al.'s (1989,
    ! S18.1) formulation.

    call eval_outer_coeffs_unno(this%bc, x_o, V_g, As, c_1)

    associate(omega_c => omega)

      lambda = outer_wavenumber(V_g, As, c_1, omega_c, 0)
      
      b_11 = V_g - 3._WP
      b_12 = -V_g

      B_o(1,1) = (lambda - b_11)/b_12
      B_o(1,2) = -1._WP

    end associate

    ! Finish

    return

  end function outer_bound_unno

!****

  function outer_bound_jcd (this, x_o, omega) result (B_o)

    class(rad_bound_t), intent(in) :: this
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

    associate(omega_c => omega)

      lambda = outer_wavenumber(V_g, As, c_1, omega, 0)

      b_11 = V_g - 3._WP
      b_12 = -V_g
      
      B_o(1,1) = (lambda - b_11)/b_12
      B_o(1,2) = -1._WP

    end associate

    ! Finish

    return

  end function outer_bound_jcd

end module gyre_rad_bound
