! Module   : gyre_bound_ad_jcd
! Purpose  : adiabatic boundary conditions (JCD formulation)
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

module gyre_bound_ad_jcd

  ! Uses

  use core_kinds

  use gyre_bound
  use gyre_coeffs
  use gyre_jacobian
  use gyre_oscpar
  use gyre_atmos

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (bound_t) :: bound_ad_jcd_t
     private
     class(coeffs_t), pointer   :: cf => null()
     class(jacobian_t), pointer :: jc => null()
     type(oscpar_t), pointer    :: op => null()
   contains 
     private
     procedure, public :: init
     procedure, public :: inner_bound
     procedure, public :: outer_bound
  end type bound_ad_jcd_t

  ! Access specifiers

  private

  public :: bound_ad_jcd_t

  ! Procedures

contains

  subroutine init (this, cf, jc, op)

    class(bound_ad_jcd_t), intent(out)    :: this
    class(coeffs_t), intent(in), target   :: cf
    class(jacobian_t), intent(in), target :: jc
    type(oscpar_t), intent(in), target    :: op

    ! Initialize the bound

    this%cf => cf
    this%jc => jc
    this%op => op

    this%n_i = 2
    this%n_o = 2
    this%n_e = this%n_i + this%n_o

    $CHECK_BOUNDS(this%n_e,this%jc%n_e)

    ! Finish

    return
    
  end subroutine init

!****

  function inner_bound (this, x_i, omega) result (B_i)

    class(bound_ad_jcd_t), intent(in) :: this
    real(WP), intent(in)              :: x_i
    complex(WP), intent(in)           :: omega
    $if($GFORTRAN_PR_58007)
    complex(WP), allocatable          :: B_i(:,:)
    $else
    complex(WP)                       :: B_i(this%n_i,this%n_e)
    $endif

    $ASSERT(x_i == 0._WP,Boundary condition invalid for x_i /= 0)

    $if($GFORTRAN_PR_58007)
    allocate(B_i(this%n_i,this%n_e))
    $endif

    ! Set the inner boundary conditions to enforce non-diverging modes

    associate(c_1 => this%cf%c_1(x_i), l => this%op%l, &
              omega_c => this%cf%omega_c(x_i, this%op%m, omega))
                 
      B_i(1,1) = 1._WP
      B_i(1,2) = -l/(c_1*omega_c**2)
      B_i(1,3) = 0._WP
      B_i(1,4) = 0._WP
        
      B_i(2,1) = 0._WP
      B_i(2,2) = 0._WP
      B_i(2,3) = l
      B_i(2,4) = -1._WP
      
    end associate

    B_i = MATMUL(B_i, this%jc%trans_matrix(x_i, omega, .TRUE.))

    ! Finish

    return

  end function inner_bound

!****

  function outer_bound (this, x_o, omega) result (B_o)

    class(bound_ad_jcd_t), intent(in) :: this
    real(WP), intent(in)              :: x_o
    complex(WP), intent(in)           :: omega
    $if($GFORTRAN_PR_58007)
    complex(WP), allocatable          :: B_o(:,:)
    $else
    complex(WP)                       :: B_o(this%n_o,this%n_e)
    $endif

    real(WP)    :: V_g
    real(WP)    :: As
    real(WP)    :: c_1
    complex(WP) :: lambda
    complex(WP) :: b_11
    complex(WP) :: b_12

    $if($GFORTRAN_PR_58007)
    allocate(B_o(this%n_o,this%n_e))
    $endif

    ! Set the outer boundary conditions

    call eval_atmos_coeffs_jcd(this%cf, x_o, V_g, As, c_1)

    associate(l => this%op%l, omega_c => this%cf%omega_c(x_o, this%op%m, omega))

      lambda = atmos_wavenumber(V_g, As, c_1, omega_c, l)

      b_11 = V_g - 3._WP
      b_12 = l*(l+1)/(c_1*omega_c**2) - V_g

      B_o(1,1) = lambda - b_11
      B_o(1,2) = -b_12
      B_o(1,3) = b_12 + (l*(l+1)/(c_1*omega_c**2) - l - 1._WP)*b_12/(V_g + As)
      B_o(1,4) = 0._WP

      B_o(2,1) = 0._WP
      B_o(2,2) = 0._WP
      B_o(2,3) = l + 1._WP
      B_o(2,4) = 1._WP

    end associate

    B_o = MATMUL(B_o, this%jc%trans_matrix(x_o, omega, .TRUE.))

    ! Finish

    return

  end function outer_bound

end module gyre_bound_ad_jcd
