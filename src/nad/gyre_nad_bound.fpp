! Incfile  : gyre_nad_bound
! Purpose  : nonadiabatic boundary conditions
!
! Copyright 2013-2017 Rich Townsend
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
  use gyre_mode_par
  use gyre_model
  use gyre_model_util
  use gyre_nad_trans
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: REGULAR_TYPE = 1
  integer, parameter :: ZERO_TYPE = 2
  integer, parameter :: DZIEM_TYPE = 3
  integer, parameter :: UNNO_TYPE = 4
  integer, parameter :: JCD_TYPE = 5

  integer, parameter :: J_V = 1
  integer, parameter :: J_V_G = 2
  integer, parameter :: J_AS = 3
  integer, parameter :: J_U = 4
  integer, parameter :: J_C_1 = 5
  integer, parameter :: J_NABLA_AD = 6

  integer, parameter :: J_LAST = J_NABLA_AD

  ! Derived-type definitions

  type, extends (c_bound_t) :: nad_bound_t
     private
     class(model_t), pointer     :: ml => null()
     class(c_rot_t), allocatable :: rt
     type(nad_trans_t)           :: tr
     real(WP), allocatable       :: coeffs(:,:)
     real(WP)                    :: alpha_gr
     complex(WP)                 :: alpha_om
     integer                     :: type_i
     integer                     :: type_o
   contains 
     private
     procedure         :: stencil_
     procedure, public :: build_i
     procedure         :: build_regular_i_
     procedure         :: build_zero_i_
     procedure, public :: build_o
     procedure         :: build_zero_o_
     procedure         :: build_dziem_o_
     procedure         :: build_unno_o_
     procedure         :: build_jcd_o_
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

  function nad_bound_t_ (ml, pt_i, pt_o, md_p, os_p) result (bd)

    class(model_t), pointer, intent(in) :: ml
    type(point_t), intent(in)           :: pt_i
    type(point_t), intent(in)           :: pt_o
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(nad_bound_t)                   :: bd

    ! Construct the nad_bound_t

    bd%ml => ml

    allocate(bd%rt, SOURCE=c_rot_t(ml, pt_i, md_p, os_p))

    bd%tr = nad_trans_t(ml, pt_i, md_p, os_p)

    select case (os_p%inner_bound)
    case ('REGULAR')
       $ASSERT(pt_i%x == 0._WP,Boundary condition invalid for x /= 0)
       bd%type_i = REGULAR_TYPE
    case ('ZERO')
       $ASSERT(pt_i%x /= 0._WP,Boundary condition invalid for x == 0)
       bd%type_i = ZERO_TYPE
    case default
       $ABORT(Invalid inner_bound)
    end select

    select case (os_p%outer_bound)
    case ('ZERO')
       bd%type_o = ZERO_TYPE
    case ('DZIEM')
       if (ml%is_vacuum(pt_o)) then
          bd%type_o = ZERO_TYPE
       else
          bd%type_o = DZIEM_TYPE
       endif
    case ('UNNO')
       if (ml%is_vacuum(pt_o)) then
          bd%type_o = ZERO_TYPE
       else
          bd%type_o = UNNO_TYPE
       end if
    case ('JCD')
       if (ml%is_vacuum(pt_o)) then
          bd%type_o = ZERO_TYPE
       else
          bd%type_o = JCD_TYPE
       end if
    case default
       $ABORT(Invalid outer_bound)
    end select

    if (os_p%cowling_approx) then
       bd%alpha_gr = 0._WP
    else
       bd%alpha_gr = 1._WP
    endif

    select case (os_p%time_factor)
    case ('OSC')
       bd%alpha_om = 1._WP
    case ('EXP')
       bd%alpha_om = (0._WP, 1._WP)
    case default
       $ABORT(Invalid time_factor)
    end select

    call bd%stencil_(pt_i, pt_o)

    bd%n_i = 3
    bd%n_o = 3

    bd%n_e = 6

    ! Finish

    return
    
  end function nad_bound_t_

  !****

  subroutine stencil_ (this, pt_i, pt_o)

    class(nad_bound_t), intent(inout) :: this
    type(point_t), intent(in)         :: pt_i
    type(point_t), intent(in)         :: pt_o

    ! Calculate coefficients at the stencil points

    call check_model(this%ml, [I_V_2,I_U,I_C_1,I_NABLA_AD])

    allocate(this%coeffs(2,J_LAST))

    ! Inner boundary

    select case (this%type_i)
    case (REGULAR_TYPE)
       this%coeffs(1,J_C_1) = this%ml%coeff(I_C_1, pt_i)
    case (ZERO_TYPE)
    case default
       $ABORT(Invalid type_i)
    end select

    ! Outer boundary

    select case (this%type_o)
    case (ZERO_TYPE)
       this%coeffs(2,J_V) = this%ml%coeff(I_V_2, pt_o)*pt_o%x**2
    case (DZIEM_TYPE)
       this%coeffs(2,J_V) = this%ml%coeff(I_V_2, pt_o)*pt_o%x**2
       this%coeffs(2,J_C_1) = this%ml%coeff(I_C_1, pt_o)
    case (UNNO_TYPE)
       call eval_atmos_coeffs_jcd(this%ml, pt_o, this%coeffs(2,J_V_G), &
            this%coeffs(2,J_AS), this%coeffs(2,J_C_1))
    case (JCD_TYPE)
       call eval_atmos_coeffs_jcd(this%ml, pt_o, this%coeffs(2,J_V_G), &
            this%coeffs(2,J_AS), this%coeffs(2,J_C_1))
    case default
       $ABORT(Invalid type_o)
    end select

    this%coeffs(2, J_U) = this%ml%coeff(I_U, pt_o)
    this%coeffs(2,J_NABLA_AD) = this%ml%coeff(I_NABLA_AD, pt_o)

    ! Set up stencils for the rt and tr components

    call this%rt%stencil([pt_i,pt_o])
    call this%tr%stencil([pt_i,pt_o])

    ! Finish

    return

  end subroutine stencil_

  !****

  subroutine build_i (this, omega, B, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    ! Evaluate the inner boundary conditions

    select case (this%type_i)
    case (REGULAR_TYPE)
       call this%build_regular_i_(omega, B, scl)
    case (ZERO_TYPE)
       call this%build_zero_i_(omega, B, scl)
    case default
       $ABORT(Invalid type_i)
    end select

    ! Finish

    return

  end subroutine build_i

  !****

  subroutine build_regular_i_ (this, omega, B, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    complex(WP) :: l_i
    complex(WP) :: omega_c

    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(scl),this%n_i)

    ! Evaluate the inner boundary conditions (regular-enforcing)

    associate( &
         c_1 => this%coeffs(1,J_C_1), &
         alpha_gr => this%alpha_gr, &
         alpha_om => this%alpha_om)

      l_i = this%rt%l_i(omega)

      omega_c = this%rt%omega_c(1, omega)

      ! Set up the boundary conditions

      B(1,1) = c_1*alpha_om*omega_c**2
      B(1,2) = -l_i
      B(1,3) = alpha_gr*(-l_i)
      B(1,4) = alpha_gr*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP

      B(2,1) = alpha_gr*(0._WP)
      B(2,2) = alpha_gr*(0._WP)
      B(2,3) = alpha_gr*(l_i)
      B(2,4) = alpha_gr*(-1._WP) + (1._WP - alpha_gr)
      B(2,5) = alpha_gr*(0._WP)
      B(2,6) = alpha_gr*(0._WP)

      B(3,1) = 0._WP
      B(3,2) = 0._WP
      B(3,3) = alpha_gr*(0._WP)
      B(3,4) = alpha_gr*(0._WP)
      B(3,5) = 1._WP
      B(3,6) = 0._WP

      scl = 1._WP

    end associate

    ! Apply the variables transformation

    call this%tr%trans_cond(B, 1, omega)

    ! Finish

    return

  end subroutine build_regular_i_

  !****

  subroutine build_zero_i_ (this, omega, B, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_i)

    ! Evaluate the inner boundary conditions (zero
    ! displacement/gravity)

    associate( &
         alpha_gr => this%alpha_gr)

      ! Set up the boundary conditions

      B(1,1) = 1._WP
      B(1,2) = 0._WP
      B(1,3) = alpha_gr*(0._WP)
      B(1,4) = alpha_gr*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP
        
      B(2,1) = alpha_gr*(0._WP)
      B(2,2) = alpha_gr*(0._WP)
      B(2,3) = alpha_gr*(0._WP)
      B(2,4) = alpha_gr*(1._WP) + (1._WP - alpha_gr)
      B(2,5) = alpha_gr*(0._WP)
      B(2,6) = alpha_gr*(0._WP)
      
      B(3,1) = 0._WP
      B(3,2) = 0._WP
      B(3,3) = alpha_gr*(0._WP)
      B(3,4) = alpha_gr*(0._WP)
      B(3,5) = 1._WP
      B(3,6) = 0._WP

      scl = 1._WP

    end associate
      
    ! Apply the variables transformation

    call this%tr%trans_cond(B, 1, omega)

    ! Finish

    return

  end subroutine build_zero_i_

  !****

  subroutine build_o (this, omega, B, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    ! Evaluate the outer boundary conditions

    select case (this%type_o)
    case (ZERO_TYPE)
       call this%build_zero_o_(omega, B, scl)
    case (DZIEM_TYPE)
       call this%build_dziem_o_(omega, B, scl)
    case (UNNO_TYPE)
       call this%build_unno_o_(omega, B, scl)
    case (JCD_TYPE)
       call this%build_jcd_o_(omega, B, scl)
    case default
       $ABORT(Invalid type_o)
    end select

    ! Finish

    return

  end subroutine build_o
  
  !****

  subroutine build_zero_o_ (this, omega, B, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    complex(WP) :: l_e

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions (zero-pressure)

    associate( &
         V => this%coeffs(2,J_V), &
         U => this%coeffs(2,J_U), &
         nabla_ad => this%coeffs(2,J_NABLA_AD), &
         alpha_gr => this%alpha_gr)

      l_e = this%rt%l_e(2, omega)

      ! Set up the boundary conditions

      B(1,1) = 1._WP
      B(1,2) = -1._WP
      B(1,3) = alpha_gr*(0._WP)
      B(1,4) = alpha_gr*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP
      
      B(2,1) = alpha_gr*(U)
      B(2,2) = alpha_gr*(0._WP)
      B(2,3) = alpha_gr*(l_e + 1._WP) + (1._WP - alpha_gr)
      B(2,4) = alpha_gr*(1._WP)
      B(2,5) = alpha_gr*(0._WP)
      B(2,6) = alpha_gr*(0._WP)

      B(3,1) = 2._WP - 4._WP*nabla_ad*V
      B(3,2) = 4._WP*nabla_ad*V
      B(3,3) = alpha_gr*(0._WP)
      B(3,4) = alpha_gr*(0._WP)
      B(3,5) = 4._WP
      B(3,6) = -1._WP

      scl = 1._WP

    end associate

    ! Apply the variables transformation

    call this%tr%trans_cond(B, 2, omega)

    ! Finish

    return

  end subroutine build_zero_o_

  !****

  subroutine build_dziem_o_ (this, omega, B, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: omega_c

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions ([Dzi1971] formulation)

    associate( &
         V => this%coeffs(2,J_V), &
         c_1 => this%coeffs(2,J_C_1), &
         nabla_ad => this%coeffs(2,J_NABLA_AD), &
         alpha_gr => this%alpha_gr, &
         alpha_om => this%alpha_om)

      lambda = this%rt%lambda(2, omega)
      l_e = this%rt%l_e(2, omega)

      omega_c = this%rt%omega_c(2, omega)

      ! Set up the boundary conditions

      B(1,1) = 1._WP + (lambda/(c_1*alpha_om*omega_c**2) - 4._WP - c_1*alpha_om*omega_c**2)/V
      B(1,2) = -1._WP
      B(1,3) = alpha_gr*((lambda/(c_1*alpha_om*omega_c**2) - l_e - 1._WP)/V)
      B(1,4) = alpha_gr*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP

      B(2,1) = alpha_gr*(0._WP)
      B(2,2) = alpha_gr*(0._WP)
      B(2,3) = alpha_gr*(l_e + 1._WP) + (1._WP - alpha_gr)
      B(2,4) = alpha_gr*(1._WP)
      B(2,5) = alpha_gr*(0._WP)
      B(2,6) = alpha_gr*(0._WP)

      B(3,1) = 2._WP - 4._WP*nabla_ad*V
      B(3,2) = 4._WP*nabla_ad*V
      B(3,3) = alpha_gr*(0._WP)
      B(3,4) = alpha_gr*(0._WP)
      B(3,5) = 4._WP
      B(3,6) = -1._WP

      scl = 1._WP

    end associate

    ! Apply the variables transformation

    call this%tr%trans_cond(B, 2, omega)

    ! Finish

    return

  end subroutine build_dziem_o_

  !****

  subroutine build_unno_o_ (this, omega, B, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: omega_c
    complex(WP) :: beta
    complex(WP) :: b_11
    complex(WP) :: b_12
    complex(WP) :: b_13
    complex(WP) :: b_21
    complex(WP) :: b_22
    complex(WP) :: b_23
    complex(WP) :: alpha_1
    complex(WP) :: alpha_2

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions ([Unn1989] formulation)

    associate( &
         V => this%coeffs(2,J_V), &
         V_g => this%coeffs(2,J_V_G), &
         As => this%coeffs(2,J_AS), &
         c_1 => this%coeffs(2,J_C_1), &
         nabla_ad => this%coeffs(2,J_NABLA_AD), &
         alpha_gr => this%alpha_gr, &
         alpha_om => this%alpha_om)

      lambda = this%rt%lambda(2, omega)
      l_e = this%rt%l_e(2, omega)

      omega_c = this%rt%omega_c(2, omega)

      beta = atmos_beta(V_g, As, c_1, omega_c, lambda)

      b_11 = V_g - 3._WP
      b_12 = lambda/(c_1*alpha_om*omega_c**2) - V_g
      b_13 = alpha_gr*(V_g)

      b_21 = c_1*alpha_om*omega_c**2 - As
      b_22 = 1._WP + As
      b_23 = alpha_gr*(-As)

      alpha_1 = (b_12*b_23 - b_13*(b_22+l_e))/((b_11+l_e)*(b_22+l_e) - b_12*b_21)
      alpha_2 = (b_21*b_13 - b_23*(b_11+l_e))/((b_11+l_e)*(b_22+l_e) - b_12*b_21)

      ! Set up the boundary conditions

      B(1,1) = beta - b_11
      B(1,2) = -b_12
      B(1,3) = -(alpha_1*(beta - b_11) - alpha_2*b_12 + b_12)
      B(1,4) = 0._WP
      B(1,5) = 0._WP
      B(1,6) = 0._WP

      B(2,1) = alpha_gr*(0._WP)
      B(2,2) = alpha_gr*(0._WP)
      B(2,3) = alpha_gr*(l_e + 1._WP) + (1._WP - alpha_gr)
      B(2,4) = alpha_gr*(1._WP)
      B(2,5) = alpha_gr*(0._WP)
      B(2,6) = alpha_gr*(0._WP)

      B(3,1) = 2._WP - 4._WP*nabla_ad*V
      B(3,2) = 4._WP*nabla_ad*V
      B(3,3) = alpha_gr*(0._WP)
      B(3,4) = alpha_gr*(0._WP)
      B(3,5) = 4._WP
      B(3,6) = -1._WP

      scl = 1._WP

    end associate

    ! Apply the variables transformation

    call this%tr%trans_cond(B, 2, omega)

    ! Finish

    return

  end subroutine build_unno_o_

  !****

  subroutine build_jcd_o_ (this, omega, B, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: omega_c
    complex(WP) :: beta
    complex(WP) :: b_11
    complex(WP) :: b_12

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions ([Chr2008] formulation)

    associate( &
         V => this%coeffs(2,J_V), &
         V_g => this%coeffs(2,J_V_G), &
         As => this%coeffs(2,J_AS), &
         c_1 => this%coeffs(2,J_C_1), &
         nabla_ad => this%coeffs(2,J_NABLA_AD), &
         alpha_gr => this%alpha_gr, &
         alpha_om => this%alpha_om)

      lambda = this%rt%lambda(2, omega)
      l_e = this%rt%l_e(2, omega)

      omega_c = this%rt%omega_c(2, omega)

      beta = atmos_beta(V_g, As, c_1, omega_c, lambda)

      b_11 = V_g - 3._WP
      b_12 = lambda/(c_1*alpha_om*omega_c**2) - V_g

      ! Set up the boundary conditions

      B(1,1) = beta - b_11
      B(1,2) = -b_12
      B(1,3) = alpha_gr*((lambda/(c_1*alpha_om*omega_c**2) - l_e - 1._WP)*b_12/(V_g + As))
      B(1,4) = alpha_gr*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP

      B(2,1) = alpha_gr*(0._WP)
      B(2,2) = alpha_gr*(0._WP)
      B(2,3) = alpha_gr*(l_e + 1._WP) + (1._WP - alpha_gr)
      B(2,4) = alpha_gr*(1._WP)
      B(2,5) = alpha_gr*(0._WP)
      B(2,6) = alpha_gr*(0._WP)

      B(3,1) = 2._WP - 4._WP*nabla_ad*V
      B(3,2) = 4._WP*nabla_ad*V
      B(3,3) = alpha_gr*(0._WP)
      B(3,4) = alpha_gr*(0._WP)
      B(3,5) = 4._WP
      B(3,6) = -1._WP

      scl = 1._WP

    end associate
    
    ! Apply the variables transformation

    call this%tr%trans_cond(B, 2, omega)

    ! Finish

    return

  end subroutine build_jcd_o_

end module gyre_nad_bound
