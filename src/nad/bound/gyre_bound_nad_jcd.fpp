! Module   : gyre_bound_nad_jcd
! Purpose  : nonadiabatic boundary conditions (JCD formulation)
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

module gyre_bound_nad_jcd

  ! Uses

  use core_kinds

  use gyre_bound
  use gyre_model
  use gyre_jacobian
  use gyre_oscpar
  use gyre_atmos

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (bound_t) :: bound_nad_jcd_t
     private
     class(model_t), pointer        :: ml => null()
     class(jacobian_t), allocatable :: jc
     type(oscpar_t)                 :: op
   contains 
     private
     procedure, public :: inner_bound => inner_bound_
     procedure, public :: outer_bound => outer_bound_
  end type bound_nad_jcd_t

  ! Interfaces

  interface bound_nad_jcd_t
     module procedure bound_nad_jcd_t_
  end interface bound_nad_jcd_t

  ! Access specifiers

  private

  public :: bound_nad_jcd_t

  ! Procedures

contains

  function bound_nad_jcd_t_ (ml, jc, op) result (bd)

    class(model_t), pointer, intent(in) :: ml
    class(jacobian_t), intent(in)       :: jc
    type(oscpar_t), intent(in)          :: op
    type(bound_nad_jcd_t)               :: bd

    ! Construct the bound_nad_jcd_t

    bd%ml => ml
    allocate(bd%jc, SOURCE=jc)
    bd%op = op

    bd%n_i = 3
    bd%n_o = 3
    bd%n_e = bd%n_i + bd%n_o

    $CHECK_BOUNDS(bd%n_e,bd%jc%n_e)

    ! Finish

    return
    
  end function bound_nad_jcd_t_

!****

  function inner_bound_ (this, x_i, omega) result (B_i)

    class(bound_nad_jcd_t), intent(in) :: this
    real(WP), intent(in)               :: x_i
    complex(WP), intent(in)            :: omega
    $if ($GFORTRAN_PR_58007)
    complex(WP), allocatable           :: B_i(:,:)
    $else
    complex(WP)                        :: B_i(this%n_i,this%n_e)
    $endif

    $ASSERT(x_i == 0._WP,Boundary condition invalid for x_i /= 0)

    $if ($GFORTRAN_PR_58007)
    allocate(B_i(this%n_i,this%n_e))
    $endif

    ! Set the inner boundary conditions to enforce non-diverging modes

    associate(c_1 => this%ml%c_1(x_i), l => this%op%l, &
              omega_c => this%ml%omega_c(x_i, this%op%m, omega))

      B_i(1,1) = c_1*omega_c**2
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

    B_i = MATMUL(B_i, this%jc%trans_matrix(x_i, omega, .TRUE.))

    ! Finish

    return

  end function inner_bound_

!****

  function outer_bound_ (this, x_o, omega) result (B_o)

    class(bound_nad_jcd_t), intent(in) :: this
    real(WP), intent(in)               :: x_o
    complex(WP), intent(in)            :: omega
    $if ($GFORTRAN_PR_58007)
    complex(WP), allocatable           :: B_o(:,:)
    $else
    complex(WP)                        :: B_o(this%n_o,this%n_e)
    $endif

    real(WP)    :: V_g
    real(WP)    :: As
    real(WP)    :: c_1
    complex(WP) :: lambda
    complex(WP) :: b_11
    complex(WP) :: b_12

    $if ($GFORTRAN_PR_58007)
    allocate(B_o(this%n_o,this%n_e))
    $endif

    ! Set the outer boundary conditions

    call eval_atmos_coeffs_jcd(this%ml, x_o, V_g, As, c_1)

    associate(V => this%ml%V(x_o), nabla_ad => this%ml%nabla_ad(x_o), &
              l => this%op%l, omega_c => this%ml%omega_c(x_o, this%op%m, omega))

      lambda = atmos_wavenumber(V_g, As, c_1, omega_c, l)

      b_11 = V_g - 3._WP
      b_12 = l*(l+1)/(c_1*omega_c**2) - V_g

      B_o(1,1) = lambda - b_11
      B_o(1,2) = -b_12
      B_o(1,3) = b_12 + (l*(l+1)/(c_1*omega_c**2) - l - 1._WP)*b_12/(V_g + As)
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

    B_o = MATMUL(B_o, this%jc%trans_matrix(x_o, omega, .TRUE.))

    ! Finish

    return

  end function outer_bound_

end module gyre_bound_nad_jcd
