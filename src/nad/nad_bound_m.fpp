! Module  : nad_bound_m
! Purpose : nonadiabatic boundary conditions
!
! Copyright 2013-2020 Rich Townsend & The GYRE Team
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

module nad_bound_m

  ! Uses

  use kinds_m

  use atmos_m
  use bound_m
  use context_m
  use math_m
  use mode_par_m
  use model_m
  use model_util_m
  use nad_trans_m
  use osc_par_m
  use point_m
  use rot_m
  use state_m

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: REGULAR_TYPE = 1
  integer, parameter :: ZERO_R_TYPE = 2
  integer, parameter :: ZERO_H_TYPE = 3
  integer, parameter :: VACUUM_TYPE = 4
  integer, parameter :: DZIEM_TYPE = 5
  integer, parameter :: ISOTHRM_TYPE = 6
  integer, parameter :: UNNO_TYPE = 7
  integer, parameter :: JCD_TYPE = 8

  integer, parameter :: J_V = 1
  integer, parameter :: J_AS = 2
  integer, parameter :: J_U = 3
  integer, parameter :: J_C_1 = 4
  integer, parameter :: J_GAMMA_1 = 5
  integer, parameter :: J_NABLA_AD = 6
  integer, parameter :: J_C_THN = 7
  integer, parameter :: J_UPS_T = 8

  integer, parameter :: J_LAST = J_UPS_T

  ! Derived-type definitions

  type, extends (c_bound_t) :: nad_bound_t
     private
     type(context_t), pointer  :: cx => null()
     type(point_t)             :: pt(2)
     type(nad_trans_t)         :: tr
     real(WP), allocatable     :: coeff(:,:)
     real(WP)                  :: alpha_grv
     real(WP)                  :: alpha_rht
     complex(WP)               :: alpha_omg
     integer                   :: type_i
     integer                   :: type_o
     character(:), allocatable :: branch_o
   contains 
     private
     procedure         :: stencil_
     procedure, public :: build_i
     procedure         :: build_regular_i_
     procedure         :: build_zero_r_i_
     procedure         :: build_zero_h_i_
     procedure, public :: build_o
     procedure         :: build_vacuum_o_
     procedure         :: build_dziem_o_
     procedure         :: build_decomp_o_
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

  function nad_bound_t_ (cx, md_p, os_p) result (bd)

    type(context_t), pointer, intent(in) :: cx
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(nad_bound_t)                    :: bd

    class(model_t), pointer :: ml
    type(point_t)           :: pt_i
    type(point_t)           :: pt_o

    ! Construct the nad_bound_t

    bd%cx => cx

    bd%tr = nad_trans_t(cx, md_p, os_p)

    ml => cx%model()
    pt_i = cx%point_i()
    pt_o = cx%point_o()

    select case (os_p%inner_bound)
    case ('REGULAR')
       $ASSERT(pt_i%x == 0._WP,Boundary condition invalid for x /= 0)
       bd%type_i = REGULAR_TYPE
    case ('ZERO_R')
       $ASSERT(pt_i%x /= 0._WP,Boundary condition invalid for x == 0)
       bd%type_i = ZERO_R_TYPE
    case ('ZERO_H')
       $ASSERT(pt_i%x /= 0._WP,Boundary condition invalid for x == 0)
       bd%type_i = ZERO_H_TYPE
    case default
       $ABORT(Invalid inner_bound)
    end select

    select case (os_p%outer_bound)
    case ('VACUUM')
       bd%type_o = VACUUM_TYPE
    case ('DZIEM')
       if (ml%is_vacuum(pt_o)) then
          bd%type_o = VACUUM_TYPE
       else
          bd%type_o = DZIEM_TYPE
       endif
    case ('ISOTHERMAL')
       if (ml%is_vacuum(pt_o)) then
          bd%type_o = VACUUM_TYPE
       else
          bd%type_o = ISOTHRM_TYPE
       end if
    case ('UNNO')
       if (ml%is_vacuum(pt_o)) then
          bd%type_o = VACUUM_TYPE
       else
          bd%type_o = UNNO_TYPE
       end if
    case ('JCD')
       if (ml%is_vacuum(pt_o)) then
          bd%type_o = VACUUM_TYPE
       else
          bd%type_o = JCD_TYPE
       end if
    case default
       $ABORT(Invalid outer_bound)
    end select

    bd%alpha_grv = os_p%alpha_grv
    bd%alpha_rht = os_p%alpha_rht
    
    select case (os_p%time_factor)
    case ('OSC')
       bd%alpha_omg = 1._WP
    case ('EXP')
       bd%alpha_omg = (0._WP, 1._WP)
    case default
       $ABORT(Invalid time_factor)
    end select

    bd%branch_o = os_p%outer_bound_branch

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

    class(model_t), pointer :: ml

    ! Calculate coefficients at the stencil points

    ml => this%cx%model()

    call check_model(ml, [I_V_2,I_U,I_C_1,I_GAMMA_1,I_NABLA_AD,I_C_THN])

    allocate(this%coeff(2,J_LAST))

    ! Inner boundary

    select case (this%type_i)
    case (REGULAR_TYPE)
       this%coeff(1,J_C_1) = ml%coeff(I_C_1, pt_i)
    case (ZERO_R_TYPE)
    case (ZERO_H_TYPE)
    case default
       $ABORT(Invalid type_i)
    end select

    ! Outer boundary

    select case (this%type_o)
    case (VACUUM_TYPE)
       this%coeff(2,J_V) = ml%coeff(I_V_2, pt_o)*pt_o%x**2
       this%coeff(2,J_U) = ml%coeff(I_U, pt_o)
    case (DZIEM_TYPE)
       this%coeff(2,J_V) = ml%coeff(I_V_2, pt_o)*pt_o%x**2
       this%coeff(2,J_C_1) = ml%coeff(I_C_1, pt_o)
    case (ISOTHRM_TYPE)
       call eval_atmos_coeffs_isothrm(ml, pt_o, this%coeff(2,J_V), &
            this%coeff(2,J_AS), this%coeff(2,J_C_1), this%coeff(2,J_GAMMA_1))
    case (UNNO_TYPE)
       call eval_atmos_coeffs_unno(ml, pt_o, this%coeff(2,J_V), &
            this%coeff(2,J_AS), this%coeff(2,J_C_1), this%coeff(2,J_GAMMA_1))
    case (JCD_TYPE)
       call eval_atmos_coeffs_isothrm(ml, pt_o, this%coeff(2,J_V), &
            this%coeff(2,J_AS), this%coeff(2,J_C_1), this%coeff(2,J_GAMMA_1))
    case default
       $ABORT(Invalid type_o)
    end select

    this%coeff(2,J_NABLA_AD) = ml%coeff(I_NABLA_AD, pt_o)
    this%coeff(2,J_C_THN) = ml%coeff(I_C_THN, pt_o)
    
    ! Set up stencil for the tr component

    call this%tr%stencil([pt_i,pt_o])

    ! Store the stencil points for on-the-fly evaluations

    this%pt = [pt_i,pt_o]

    ! Finish

    return

  end subroutine stencil_

  !****

  subroutine build_i (this, st, B, scl)

    class(nad_bound_t), intent(in) :: this
    class(c_state_t), intent(in)   :: st
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    ! Evaluate the inner boundary conditions

    select case (this%type_i)
    case (REGULAR_TYPE)
       call this%build_regular_i_(st, B, scl)
    case (ZERO_R_TYPE)
       call this%build_zero_r_i_(st, B, scl)
    case (ZERO_H_TYPE)
       call this%build_zero_h_i_(st, B, scl)
    case default
       $ABORT(Invalid type_i)
    end select

    ! Apply the variables transformation

    call this%tr%trans_cond(B, 1, st)

    ! Finish

    return

  end subroutine build_i

  !****

  subroutine build_regular_i_ (this, st, B, scl)

    class(nad_bound_t), intent(in) :: this
    class(c_state_t), intent(in)   :: st
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    real(WP)    :: Omega_rot
    complex(WP) :: omega_c
    complex(WP) :: l_i

    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(scl),this%n_i)

    ! Evaluate the inner boundary conditions (regular-enforcing)

    associate( &
         c_1 => this%coeff(1,J_C_1), &
         pt => this%pt(1), &
         alpha_grv => this%alpha_grv, &
         alpha_omg => this%alpha_omg)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, st)

      l_i = this%cx%l_e(Omega_rot, st)

      ! Set up the boundary conditions

      B(1,1) = c_1*alpha_omg*omega_c**2
      B(1,2) = -l_i
      B(1,3) = alpha_grv*(-l_i)
      B(1,4) = alpha_grv*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP

      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(l_i)
      B(2,4) = alpha_grv*(-1._WP) + (1._WP - alpha_grv)
      B(2,5) = alpha_grv*(0._WP)
      B(2,6) = alpha_grv*(0._WP)

      B(3,1) = 0._WP
      B(3,2) = 0._WP
      B(3,3) = alpha_grv*(0._WP)
      B(3,4) = alpha_grv*(0._WP)
      B(3,5) = 1._WP
      B(3,6) = 0._WP

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_regular_i_

  !****

  subroutine build_zero_r_i_ (this, st, B, scl)

    class(nad_bound_t), intent(in) :: this
    class(c_state_t), intent(in)   :: st
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_i)

    ! Evaluate the inner boundary conditions (zero
    ! radial displacement/gravity)

    associate( &
         alpha_grv => this%alpha_grv)

      ! Set up the boundary conditions

      B(1,1) = 1._WP
      B(1,2) = 0._WP
      B(1,3) = alpha_grv*(0._WP)
      B(1,4) = alpha_grv*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP
        
      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(0._WP)
      B(2,4) = alpha_grv*(1._WP) + (1._WP - alpha_grv)
      B(2,5) = alpha_grv*(0._WP)
      B(2,6) = alpha_grv*(0._WP)
      
      B(3,1) = 0._WP
      B(3,2) = 0._WP
      B(3,3) = alpha_grv*(0._WP)
      B(3,4) = alpha_grv*(0._WP)
      B(3,5) = 1._WP
      B(3,6) = 0._WP

      scl = 1._WP

    end associate
      
    ! Finish

    return

  end subroutine build_zero_r_i_

  !****

  subroutine build_zero_h_i_ (this, st, B, scl)

    class(nad_bound_t), intent(in) :: this
    class(c_state_t), intent(in)   :: st
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_i)

    ! Evaluate the inner boundary conditions (zero
    ! horizontal displacement/gravity)

    associate( &
         alpha_grv => this%alpha_grv)

      ! Set up the boundary conditions

      B(1,1) = 0._WP
      B(1,2) = 1._WP
      B(1,3) = alpha_grv*(1._WP)
      B(1,4) = alpha_grv*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP
        
      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(0._WP)
      B(2,4) = alpha_grv*(1._WP) + (1._WP - alpha_grv)
      B(2,5) = alpha_grv*(0._WP)
      B(2,6) = alpha_grv*(0._WP)
      
      B(3,1) = 0._WP
      B(3,2) = 0._WP
      B(3,3) = alpha_grv*(0._WP)
      B(3,4) = alpha_grv*(0._WP)
      B(3,5) = 1._WP
      B(3,6) = 0._WP

      scl = 1._WP

    end associate
      
    ! Finish

    return

  end subroutine build_zero_h_i_

  !****

  subroutine build_o (this, st, B, scl)

    class(nad_bound_t), intent(in) :: this
    class(c_state_t), intent(in)   :: st
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    ! Evaluate the outer boundary conditions

    select case (this%type_o)
    case (VACUUM_TYPE)
       call this%build_vacuum_o_(st, B, scl)
    case (DZIEM_TYPE)
       call this%build_dziem_o_(st, B, scl)
    case (ISOTHRM_TYPE)
       call this%build_decomp_o_(st, B, scl)
    case (UNNO_TYPE)
       call this%build_decomp_o_(st, B, scl)
    case (JCD_TYPE)
       call this%build_jcd_o_(st, B, scl)
    case default
       $ABORT(Invalid type_o)
    end select

    ! Apply the variables transformation

    call this%tr%trans_cond(B, 2, st)

    ! Finish

    return

  end subroutine build_o
  
  !****

  subroutine build_vacuum_o_ (this, st, B, scl)

    class(nad_bound_t), intent(in) :: this
    class(c_state_t), intent(in)   :: st
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    real(WP)    :: Omega_rot
    complex(WP) :: omega_c
    complex(WP) :: i_omega_c
    complex(WP) :: l_e
    complex(WP) :: f_rht

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions (vacuum)

    associate( &
         V => this%coeff(2,J_V), &
         U => this%coeff(2,J_U), &
         nabla_ad => this%coeff(2,J_NABLA_AD), &
         c_thn => this%coeff(2,J_C_THN), &
         pt => this%pt(2), &
         alpha_grv => this%alpha_grv, &
         alpha_rht => this%alpha_rht, &
         alpha_omg => this%alpha_omg)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, st)
      i_omega_c = (0._WP,1._WP)*sqrt(CMPLX(alpha_omg, KIND=WP))*omega_c

      l_e = this%cx%l_e(Omega_rot, st)

      f_rht = 1._WP - 0.25_WP*alpha_rht*i_omega_c*c_thn

      ! Set up the boundary conditions

      B(1,1) = 1._WP
      B(1,2) = -1._WP
      B(1,3) = alpha_grv*(0._WP)
      B(1,4) = alpha_grv*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP
      
      B(2,1) = alpha_grv*(U)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(l_e + 1._WP) + (1._WP - alpha_grv)
      B(2,4) = alpha_grv*(1._WP)
      B(2,5) = alpha_grv*(0._WP)
      B(2,6) = alpha_grv*(0._WP)

      B(3,1) = 2._WP - 4._WP*nabla_ad*V
      B(3,2) = 4._WP*nabla_ad*V
      B(3,3) = alpha_grv*(0._WP)
      B(3,4) = alpha_grv*(0._WP)
      B(3,5) = 4._WP*f_rht
      B(3,6) = -1._WP

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_vacuum_o_

  !****

  subroutine build_dziem_o_ (this, st, B, scl)

    class(nad_bound_t), intent(in) :: this
    class(c_state_t), intent(in)   :: st
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    real(WP)    :: Omega_rot
    complex(WP) :: omega_c
    complex(WP) :: i_omega_c
    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: f_rht
    
    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions ([Dzi1971] formulation)

    associate( &
         V => this%coeff(2,J_V), &
         c_1 => this%coeff(2,J_C_1), &
         nabla_ad => this%coeff(2,J_NABLA_AD), &
         c_thn => this%coeff(2,J_C_THN), &
         pt => this%pt(2), &
         alpha_grv => this%alpha_grv, &
         alpha_rht => this%alpha_rht, &
         alpha_omg => this%alpha_omg)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, st)
      i_omega_c = (0._WP,1._WP)*sqrt(CMPLX(alpha_omg, KIND=WP))*omega_c

      lambda = this%cx%lambda(Omega_rot, st)
      l_e = this%cx%l_e(Omega_rot, st)

      f_rht = 1._WP - 0.25_WP*alpha_rht*i_omega_c*c_thn

      ! Set up the boundary conditions

      B(1,1) = 1._WP + (lambda/(c_1*alpha_omg*omega_c**2) - 4._WP - c_1*alpha_omg*omega_c**2)/V
      B(1,2) = -1._WP
      B(1,3) = alpha_grv*((lambda/(c_1*alpha_omg*omega_c**2) - l_e - 1._WP)/V)
      B(1,4) = alpha_grv*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP

      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(l_e + 1._WP) + (1._WP - alpha_grv)
      B(2,4) = alpha_grv*(1._WP)
      B(2,5) = alpha_grv*(0._WP)
      B(2,6) = alpha_grv*(0._WP)

      B(3,1) = 2._WP - 4._WP*nabla_ad*V
      B(3,2) = 4._WP*nabla_ad*V
      B(3,3) = alpha_grv*(0._WP)
      B(3,4) = alpha_grv*(0._WP)
      B(3,5) = 4._WP*f_rht
      B(3,6) = -1._WP

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_dziem_o_

  !****

  subroutine build_decomp_o_ (this, st, B, scl)

    class(nad_bound_t), intent(in) :: this
    class(c_state_t), intent(in)   :: st
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    real(WP)    :: Omega_rot
    complex(WP) :: omega_c
    complex(WP) :: i_omega_c
    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: f_rht
    complex(WP) :: a_11
    complex(WP) :: a_12
    complex(WP) :: a_13
    complex(WP) :: a_14
    complex(WP) :: a_21
    complex(WP) :: a_22
    complex(WP) :: a_23
    complex(WP) :: a_24
    complex(WP) :: chi
    complex(WP) :: G_1
    complex(WP) :: G_2

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions, based on a local
    ! eigendecomposition

    associate( &
         V => this%coeff(2,J_V), &
         As => this%coeff(2,J_AS), &
         c_1 => this%coeff(2,J_C_1), &
         Gamma_1 => this%coeff(2,J_GAMMA_1), &
         nabla_ad => this%coeff(2,J_NABLA_AD), &
         c_thn => this%coeff(2,J_C_THN), &
         pt => this%pt(2), &
         alpha_grv => this%alpha_grv, &
         alpha_rht => this%alpha_rht, &
         alpha_omg => this%alpha_omg, &
         branch => this%branch_o)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, st)
      i_omega_c = (0._WP,1._WP)*sqrt(CMPLX(alpha_omg, KIND=WP))*omega_c

      lambda = this%cx%lambda(Omega_rot, st)
      l_e = this%cx%l_e(Omega_rot, st)

      f_rht = 1._WP - 0.25_WP*alpha_rht*i_omega_c*c_thn
      
      ! Evaluate selected elements of the Jacobian matrix

      a_11 = V/Gamma_1 - 3._WP
      a_12 = lambda/(c_1*alpha_omg*omega_c**2) - V/Gamma_1
      a_13 = alpha_grv*(lambda/(c_1*alpha_omg*omega_c**2))
      a_14 = alpha_grv*(0._WP)

      a_21 = c_1*alpha_omg*omega_c**2 - As
      a_22 = As + 1._WP
      a_23 = alpha_grv*(0._WP)
      a_24 = alpha_grv*(-1._WP)

      ! Evaluate the eigenvalue for the wave we want to keep

      chi = atmos_chi(V, As, c_1, Gamma_1, omega_c, lambda, branch)

      ! Evaluate selected elements of the left eigenvectors

      G_1 = (l_e*H_(-l_e, chi) + (l_e+1._WP)*H_(l_e+1._WP, chi))/(2._WP*l_e + 1._WP)
      G_2 = (H_(-l_e, chi) - H_(l_e+1._WP, chi))/(2._WP*l_e + 1._WP)

      ! Set up the boundary conditions

      B(1,1) = -(chi - a_11)
      B(1,2) = a_12
      B(1,3) = -alpha_grv*G_1
      B(1,4) = alpha_grv*G_2
      B(1,5) = 0._WP
      B(1,6) = 0._WP

      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(l_e + 1._WP) + (1._WP - alpha_grv)
      B(2,4) = alpha_grv*(1._WP)
      B(2,5) = alpha_grv*(0._WP)
      B(2,6) = alpha_grv*(0._WP)

      B(3,1) = 2._WP - 4._WP*nabla_ad*V
      B(3,2) = 4._WP*nabla_ad*V
      B(3,3) = alpha_grv*(0._WP)
      B(3,4) = alpha_grv*(0._WP)
      B(3,5) = 4._WP*f_rht
      B(3,6) = -1._WP

      scl = 1._WP

    end associate

    ! Finish

    return

  contains

    function H_ (l, chi) result (H)

      complex(WP), intent(in) :: l
      complex(WP), intent(in) :: chi
      complex(WP)             :: H

      ! Evaluate the H function defined in leaky boundary condition
      ! notes

      H = ((a_11 - chi)*(a_13 - a_14*(1._WP-l)) + a_12*(a_23 - a_24*(1._WP-l))) / &
          (chi + l - a_11 - a_12)

      ! Finish

      return

    end function H_

  end subroutine build_decomp_o_

  !****

  subroutine build_jcd_o_ (this, st, B, scl)

    class(nad_bound_t), intent(in) :: this
    class(c_state_t), intent(in)   :: st
    complex(WP), intent(out)       :: B(:,:)
    complex(WP), intent(out)       :: scl(:)

    real(WP)    :: Omega_rot
    complex(WP) :: omega_c
    complex(WP) :: i_omega_c
    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: f_rht
    complex(WP) :: chi
    complex(WP) :: b_11
    complex(WP) :: b_12

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions ([Chr2008] formulation)

    associate( &
         V => this%coeff(2,J_V), &
         As => this%coeff(2,J_AS), &
         c_1 => this%coeff(2,J_C_1), &
         Gamma_1 => this%coeff(2,J_GAMMA_1), &
         nabla_ad => this%coeff(2,J_NABLA_AD), &
         c_thn => this%coeff(2,J_C_THN), &
         pt => this%pt(2), &
         alpha_grv => this%alpha_grv, &
         alpha_rht => this%alpha_rht, &
         alpha_omg => this%alpha_omg, &
         branch => this%branch_o)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, st)
      i_omega_c = (0._WP,1._WP)*sqrt(CMPLX(alpha_omg, KIND=WP))*omega_c

      lambda = this%cx%lambda(Omega_rot, st)
      l_e = this%cx%l_e(Omega_rot, st)

      f_rht = 1._WP - 0.25_WP*alpha_rht*i_omega_c*c_thn

      chi = atmos_chi(V, As, c_1, Gamma_1, omega_c, lambda, branch)

      b_11 = V/Gamma_1 - 3._WP
      b_12 = lambda/(c_1*alpha_omg*omega_c**2) - V/Gamma_1

      ! Set up the boundary conditions

      B(1,1) = chi - b_11
      B(1,2) = -b_12
      B(1,3) = alpha_grv*((lambda/(c_1*alpha_omg*omega_c**2) - l_e - 1._WP)*b_12/(V/Gamma_1 + As))
      B(1,4) = alpha_grv*(0._WP)
      B(1,5) = 0._WP
      B(1,6) = 0._WP

      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(l_e + 1._WP) + (1._WP - alpha_grv)
      B(2,4) = alpha_grv*(1._WP)
      B(2,5) = alpha_grv*(0._WP)
      B(2,6) = alpha_grv*(0._WP)

      B(3,1) = 2._WP - 4._WP*nabla_ad*V
      B(3,2) = 4._WP*nabla_ad*V
      B(3,3) = alpha_grv*(0._WP)
      B(3,4) = alpha_grv*(0._WP)
      B(3,5) = 4._WP*f_rht
      B(3,6) = -1._WP

      scl = 1._WP

    end associate
    
    ! Finish

    return

  end subroutine build_jcd_o_

end module nad_bound_m
