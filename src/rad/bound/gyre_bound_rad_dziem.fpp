! Module   : gyre_bound_rad_dziem
! Purpose  : radial adiabatic boundary conditions (Dziembowski formulation)
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

module gyre_bound_rad_dziem

  ! Uses

  use core_kinds

  use gyre_bound
  use gyre_coeffs
  use gyre_jacobian
  use gyre_oscpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (bound_t) :: bound_rad_dziem_t
     private
     class(coeffs_t), pointer   :: cf => null()
     class(jacobian_t), pointer :: jc => null()
     type(oscpar_t), pointer    :: op => null()
   contains 
     private
     procedure, public :: init
     procedure, public :: inner_bound
     procedure, public :: outer_bound
  end type bound_rad_dziem_t
  
  ! Access specifiers

  private

  public :: bound_rad_dziem_t

  ! Procedures

contains

  subroutine init (this, cf, jc, op)

    class(bound_rad_dziem_t), intent(out) :: this
    class(coeffs_t), intent(in), target   :: cf
    class(jacobian_t), intent(in), target :: jc
    type(oscpar_t), intent(in), target    :: op

    ! Initialize the bound

    this%cf => cf
    this%op => op

    this%n_i = 1
    this%n_o = 1
    this%n_e = this%n_i + this%n_o

    $CHECK_BOUNDS(this%n_e,this%jc%n_e)

    ! Finish

    return
    
  end subroutine init

!****

  function inner_bound (this, x_i, omega) result (B_i)

    class(bound_rad_dziem_t), intent(in) :: this
    real(WP), intent(in)                 :: x_i
    complex(WP), intent(in)              :: omega
    complex(WP)                          :: B_i(this%n_i,this%n_e)

    $ASSERT(x_i == 0._WP,Boundary condition invalid for x_i /= 0)

    ! Set the inner boundary conditions to enforce non-diverging modes

    B_i(1,1) = 1._WP
    B_i(1,2) = 0._WP

    B_i = MATMUL(B_i, this%jc%trans_matrix(x_i, omega, .TRUE.))

    ! Finish

    return

  end function inner_bound

!****

  function outer_bound (this, x_o, omega) result (B_o)

    class(bound_rad_dziem_t), intent(in) :: this
    real(WP), intent(in)                 :: x_o
    complex(WP), intent(in)              :: omega
    complex(WP)                          :: B_o(this%n_o,this%n_e)

    ! Set the outer boundary conditions

    associate(V => this%cf%V(x_o), c_1 => this%cf%c_1(x_o), &
              omega_c => omega)
        
      B_o(1,1) = 1 - (4._WP + c_1*omega_c**2)/V
      B_o(1,2) = -1._WP

    end associate

    B_o = MATMUL(B_o, this%jc%trans_matrix(x_o, omega, .TRUE.))

    ! Finish

    return

  end function outer_bound

end module gyre_bound_rad_dziem
