! Incfile  : gyre_nad_bound
! Purpose  : boundary conditions (nonadiabatic)
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

module gyre_nad_bound

  ! Uses

  use core_kinds

  use gyre_atmos
  use gyre_bound
  use gyre_eqns
  use gyre_model
  use gyre_osc_par
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

  type, extends (c_bound_t) :: nad_bound_t
     private
     class(model_t), pointer      :: ml => null()
     class(c_rot_t), allocatable  :: rt
     class(c_eqns_t), allocatable :: eq
     real(WP)                     :: x_i
     real(WP)                     :: x_o
     integer                      :: type_i
     integer                      :: type_o
     logical                      :: cowling_approx
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
  end type nad_bound_t

  ! Interfaces

  interface nad_bound_t
     module procedure nad_bound_t_
  end interface nad_bound_t

  ! Access specifiers

  private

  public :: nad_bound_t

  ! Procedures

contains

  function nad_bound_t_ (ml, rt, eq, op, x_i, x_o) result (bd)

    class(model_t), pointer, intent(in) :: ml
    class(c_rot_t), intent(in)          :: rt
    class(c_eqns_t), intent(in)         :: eq
    type(osc_par_t), intent(in)         :: op
    real(WP)                            :: x_i
    real(WP)                            :: x_o
    type(nad_bound_t)                   :: bd

    ! Construct the nad_bound_t

    bd%ml => ml
    allocate(bd%rt, SOURCE=rt)
    allocate(bd%eq, SOURCE=eq)

    bd%x_i = x_i
    bd%x_o = x_o

    select case (op%inner_bound)
    case ('REGULAR')
       bd%type_i = REGULAR_TYPE_I
    case ('ZERO')
       bd%type_i = ZERO_TYPE_I
    case default
       $ABORT(Invalid inner_bound)
    end select

    select case (op%outer_bound)
    case ('ZERO')
       bd%type_o = ZERO_TYPE_O
    case ('DZIEM')
       bd%type_o = DZIEM_TYPE_O
    case ('UNNO')
       bd%type_o = UNNO_TYPE_O
    case ('JCD')
       bd%type_o = JCD_TYPE_O
    case default
       $ABORT(Invalid outer_bound)
    end select

    bd%cowling_approx = op%cowling_approx

    bd%n_i = 3
    bd%n_o = 3
    bd%n_e = bd%n_i + bd%n_o

    ! Finish

    return
    
  end function nad_bound_t_

!****

  function B_i_ (this, omega) result (B_i)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_i(this%n_i,this%n_e)

    ! Evaluate the inner boundary conditions

    select case (this%type_i)
    case (REGULAR_TYPE_I)
       B_i = this%B_i_regular_(omega)
    case (ZERO_TYPE_I)
       B_i = this%B_i_zero_(omega)
    case default
       $ABORT(Invalid type_i)
    end select

    ! Transform to the variables used in the differential equations

    B_i = MATMUL(B_i, this%eq%T(this%x_i, omega, .TRUE.))

    ! Finish

    return

  end function B_i_

!****

  function B_i_regular_ (this, omega) result (B_i)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_i(this%n_i,this%n_e)

    real(WP)    :: c_1
    complex(WP) :: l_e
    complex(WP) :: omega_c
    real(WP)    :: chi_cowl

    $ASSERT(this%x_i == 0._WP,Boundary condition invalid for x_i /= 0)

    ! Evaluate the inner boundary conditions (regular-enforcing)

    ! Calculate coefficients

    c_1 = this%ml%c_1(this%x_i)

    l_e = this%rt%l_e(this%x_i, omega)

    omega_c = this%rt%omega_c(this%x_i, omega)

    if (this%cowling_approx) then
       chi_cowl = 0._WP
    else
       chi_cowl = 1._WP
    endif

    ! Set up the boundary conditions

    B_i(1,1) = c_1*omega_c**2
    B_i(1,2) = -l_e
    B_i(1,3) = chi_cowl*(0._WP)
    B_i(1,4) = chi_cowl*(0._WP)
    B_i(1,5) = 0._WP
    B_i(1,6) = 0._WP

    B_i(2,1) = chi_cowl*(0._WP)
    B_i(2,2) = chi_cowl*(0._WP)
    B_i(2,3) = chi_cowl*(l_e)
    B_i(2,4) = chi_cowl*(-1._WP) + (1._WP - chi_cowl)
    B_i(2,5) = chi_cowl*(0._WP)
    B_i(2,6) = chi_cowl*(0._WP)

    B_i(3,1) = 0._WP
    B_i(3,2) = 0._WP
    B_i(3,3) = chi_cowl*(0._WP)
    B_i(3,4) = chi_cowl*(0._WP)
    B_i(3,5) = 1._WP
    B_i(3,6) = 0._WP

    ! Finish

    return

  end function B_i_regular_

!****

  function B_i_zero_ (this, omega) result (B_i)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_i(this%n_i,this%n_e)

    real(WP)    :: chi_cowl

    $ASSERT(this%x_i /= 0._WP,Boundary condition invalid for x_i == 0)

    ! Evaluate the inner boundary conditions (zero
    ! displacement/gravity/entropy)

    ! Calculate coefficients

    if (this%cowling_approx) then
       chi_cowl = 0._WP
    else
       chi_cowl = 1._WP
    endif

    ! Set up the boundary conditions

    B_i(1,1) = 1._WP
    B_i(1,2) = 0._WP
    B_i(1,3) = chi_cowl*(0._WP)
    B_i(1,4) = chi_cowl*(0._WP)
    B_i(1,5) = 0._WP
    B_i(1,6) = 0._WP
        
    B_i(2,1) = chi_cowl*(0._WP)
    B_i(2,2) = chi_cowl*(0._WP)
    B_i(2,3) = chi_cowl*(0._WP)
    B_i(2,4) = chi_cowl*(1._WP) + (1._WP - chi_cowl)
    B_i(2,5) = chi_cowl*(0._WP)
    B_i(2,6) = chi_cowl*(0._WP)
      
    B_i(3,1) = 0._WP
    B_i(3,2) = 0._WP
    B_i(3,3) = chi_cowl*(0._WP)
    B_i(3,4) = chi_cowl*(0._WP)
    B_i(3,5) = 1._WP
    B_i(3,6) = 0._WP
      
    ! Finish

    return

  end function B_i_zero_

!****

  function B_o_ (this, omega) result (B_o)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

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

    ! Transform to the variables used in the differential equations

    B_o = MATMUL(B_o, this%eq%T(this%x_o, omega, .TRUE.))
    
    ! Finish

    return

  end function B_o_

!****

  function B_o_zero_ (this, omega) result (B_o)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    real(WP)    :: V
    real(WP)    :: U
    real(WP)    :: nabla_ad
    complex(WP) :: l_e
    real(WP)    :: chi_cowl

    ! Evaluate the outer boundary conditions (zero-pressure)

    ! Calculate coefficients

    V = this%ml%V_2(this%x_o)*this%x_o**2
    U = this%ml%U(this%x_o)
    nabla_ad = this%ml%nabla_ad(this%x_o)

    l_e = this%rt%l_e(this%x_i, omega)

    if (this%cowling_approx) then
       chi_cowl = 0._WP
    else
       chi_cowl = 1._WP
    endif

    ! Set up the boundary conditions

    B_o(1,1) = 1._WP
    B_o(1,2) = -1._WP
    B_o(1,3) = chi_cowl*(1._WP)
    B_o(1,4) = chi_cowl*(0._WP)
    B_o(1,5) = 0._WP
    B_o(1,6) = 0._WP
      
    B_o(2,1) = chi_cowl*(U)
    B_o(2,2) = chi_cowl*(0._WP)
    B_o(2,3) = chi_cowl*(l_e + 1._WP) + (1._WP - chi_cowl)
    B_o(2,4) = chi_cowl*(1._WP)
    B_o(2,5) = chi_cowl*(0._WP)
    B_o(2,6) = chi_cowl*(0._WP)

    B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
    B_o(3,2) = 4._WP*nabla_ad*V
    B_o(3,3) = chi_cowl*(-4._WP*nabla_ad*V)
    B_o(3,4) = chi_cowl*(0._WP)
    B_o(3,5) = 4._WP
    B_o(3,6) = -1._WP

    ! Finish

    return

  end function B_o_zero_

!****

  function B_o_dziem_ (this, omega) result (B_o)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    real(WP)    :: V
    real(WP)    :: nabla_ad
    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: omega_c
    real(WP)    :: chi_cowl

    ! Evaluate the outer boundary conditions ([Dzi1971] formulation)

    ! Calculate coefficients

    V = this%ml%V_2(this%x_o)*this%x_o**2
    nabla_ad = this%ml%nabla_ad(this%x_o)

    lambda = this%rt%lambda(this%x_o, omega)
    l_e = this%rt%l_e(this%x_o, omega)
    
    omega_c = this%rt%omega_c(this%x_o, omega)

    if (this%cowling_approx) then
       chi_cowl = 0._WP
    else
       chi_cowl = 1._WP
    endif

    ! Set up the boundary conditions

    B_o(1,1) = 1 + (lambda/omega_c**2 - 4 - omega_c**2)/V
    B_o(1,2) = -1._WP
    B_o(1,3) = chi_cowl*(1 + (lambda/omega_c**2 - l_e - 1)/V)
    B_o(1,4) = chi_cowl*(0._WP)
    B_o(1,5) = 0._WP
    B_o(1,6) = 0._WP
     
    B_o(2,1) = chi_cowl*(0._WP)
    B_o(2,2) = chi_cowl*(0._WP)
    B_o(2,3) = chi_cowl*(l_e + 1._WP) + (1._WP - chi_cowl)
    B_o(2,4) = chi_cowl*(1._WP)
    B_o(2,5) = chi_cowl*(0._WP)
    B_o(2,6) = chi_cowl*(0._WP)

    B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
    B_o(3,2) = 4._WP*nabla_ad*V
    B_o(3,3) = chi_cowl*(-4._WP*nabla_ad*V)
    B_o(3,4) = chi_cowl*(0._WP)
    B_o(3,5) = 4._WP
    B_o(3,6) = -1._WP

    ! Finish

    return

  end function B_o_dziem_

!****

  function B_o_unno_ (this, omega) result (B_o)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    real(WP)    :: V_g
    real(WP)    :: As
    real(WP)    :: c_1
    real(WP)    :: V
    real(WP)    :: nabla_ad
    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: omega_c
    complex(WP) :: beta
    real(WP)    :: chi_cowl
    complex(WP) :: b_11
    complex(WP) :: b_12
    complex(WP) :: b_13
    complex(WP) :: b_21
    complex(WP) :: b_22
    complex(WP) :: b_23
    complex(WP) :: alpha_1
    complex(WP) :: alpha_2

    ! Evaluate the outer boundary conditions ([Unn1989] formulation)

    ! Calculate coefficients

    call eval_atmos_coeffs_unno(this%ml, this%x_o, V_g, As, c_1)

    V = this%ml%V_2(this%x_o)*this%x_o**2
    nabla_ad = this%ml%nabla_ad(this%x_o)

    lambda = this%rt%lambda(this%x_o, omega)
    l_e = this%rt%l_e(this%x_o, omega)

    omega_c = this%rt%omega_c(this%x_o, omega)

    beta = atmos_beta(V_g, As, c_1, omega_c, lambda)

    if (this%cowling_approx) then
       chi_cowl = 0._WP
    else
       chi_cowl = 1._WP
    endif

    b_11 = V_g - 3._WP
    b_12 = lambda/(c_1*omega_c**2) - V_g
    b_13 = chi_cowl*(V_g)

    b_21 = c_1*omega_c**2 - As
    b_22 = 1._WP + As
    b_23 = chi_cowl*(-As)
    
    alpha_1 = (b_12*b_23 - b_13*(b_22+l_e))/((b_11+l_e)*(b_22+l_e) - b_12*b_21)
    alpha_2 = (b_21*b_13 - b_23*(b_11+l_e))/((b_11+l_e)*(b_22+l_e) - b_12*b_21)

    ! Set up the boundary conditions

    B_o(1,1) = beta - b_11
    B_o(1,2) = -b_12
    B_o(1,3) = -(alpha_1*(beta - b_11) - alpha_2*b_12)
    B_o(1,4) = 0._WP
    B_o(1,5) = 0._WP
    B_o(1,6) = 0._WP

    B_o(2,1) = chi_cowl*(0._WP)
    B_o(2,2) = chi_cowl*(0._WP)
    B_o(2,3) = chi_cowl*(l_e + 1._WP) + (1._WP - chi_cowl)
    B_o(2,4) = chi_cowl*(1._WP)
    B_o(2,5) = chi_cowl*(0._WP)
    B_o(2,6) = chi_cowl*(0._WP)
    
    B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
    B_o(3,2) = 4._WP*nabla_ad*V
    B_o(3,3) = chi_cowl*(-4._WP*nabla_ad*V)
    B_o(3,4) = chi_cowl*(0._WP)
    B_o(3,5) = 4._WP
    B_o(3,6) = -1._WP

    ! Finish

    return

  end function B_o_unno_

!****

  function B_o_jcd_ (this, omega) result (B_o)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: B_o(this%n_o,this%n_e)

    real(WP)    :: V_g
    real(WP)    :: As
    real(WP)    :: c_1
    real(WP)    :: V
    real(WP)    :: nabla_ad
    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: omega_c
    complex(WP) :: beta
    real(WP)    :: chi_cowl
    complex(WP) :: b_11
    complex(WP) :: b_12

    ! Evaluate the outer boundary conditions ([Chr2008] formulation)

    ! Calculate coefficients

    call eval_atmos_coeffs_jcd(this%ml, this%x_o, V_g, As, c_1)

    V = this%ml%V_2(this%x_o)*this%x_o**2
    nabla_ad = this%ml%nabla_ad(this%x_o)
    
    lambda = this%rt%lambda(this%x_o, omega)
    l_e = this%rt%l_e(this%x_o, omega)

    omega_c = this%rt%omega_c(this%x_o, omega)

    beta = atmos_beta(V_g, As, c_1, omega_c, lambda)

    if (this%cowling_approx) then
       chi_cowl = 0._WP
    else
       chi_cowl = 1._WP
    endif

    b_11 = V_g - 3._WP
    b_12 = lambda/(c_1*omega_c**2) - V_g

    ! Set up the boundary conditions

    B_o(1,1) = beta - b_11
    B_o(1,2) = -b_12
    B_o(1,3) = chi_cowl*(b_12 + (lambda/(c_1*omega_c**2) - l_e - 1._WP)*b_12/(V_g + As))
    B_o(1,4) = chi_cowl*(0._WP)
    B_o(1,5) = 0._WP
    B_o(1,6) = 0._WP
    
    B_o(2,1) = chi_cowl*(0._WP)
    B_o(2,2) = chi_cowl*(0._WP)
    B_o(2,3) = chi_cowl*(l_e + 1._WP) + (1._WP - chi_cowl)
    B_o(2,4) = chi_cowl*(1._WP)
    B_o(2,5) = chi_cowl*(0._WP)
    B_o(2,6) = chi_cowl*(0._WP)
    
    B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
    B_o(3,2) = 4._WP*nabla_ad*V
    B_o(3,3) = chi_cowl*(-4._WP*nabla_ad*V)
    B_o(3,4) = chi_cowl*(0._WP)
    B_o(3,5) = 4._WP
    B_o(3,6) = -1._WP

    ! Finish

    return

  end function B_o_jcd_

end module gyre_nad_bound
