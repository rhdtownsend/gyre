! Incfile  : gyre_ad_bound
! Purpose  : adiabatic boundary conditions
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

module gyre_ad_bound

  ! Uses

  use core_kinds

  use gyre_ad_trans
  use gyre_atmos
  use gyre_bound
  use gyre_context
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_state

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
  integer, parameter :: GAMMA_TYPE = 9

  integer, parameter :: J_V = 1
  integer, parameter :: J_AS = 2
  integer, parameter :: J_U = 3
  integer, parameter :: J_C_1 = 4
  integer, parameter :: J_GAMMA_1 = 5

  integer, parameter :: J_LAST = J_GAMMA_1

  ! Derived-type definitions

  type, extends (r_bound_t) :: ad_bound_t
     private
     type(context_t), pointer :: cx => null()
     type(point_t)            :: pt(2)
     type(ad_trans_t)         :: tr
     real(WP), allocatable    :: coeff(:,:)
     real(WP)                 :: alpha_grv
     real(WP)                 :: alpha_omg
     integer                  :: type_i
     integer                  :: type_o
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
     procedure         :: build_gamma_o_
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

  function ad_bound_t_ (cx, md_p, os_p) result (bd)

    type(context_t), pointer, intent(in) :: cx
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(ad_bound_t)                     :: bd

    class(model_t), pointer :: ml
    type(point_t)           :: pt_i
    type(point_t)           :: pt_o

    ! Construct the ad_bound_t

    bd%cx => cx
    
    bd%tr = ad_trans_t(cx, md_p, os_p)

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
       end if
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
    case ('GAMMA')
        bd%type_o = GAMMA_TYPE
    case default
       $ABORT(Invalid outer_bound)
    end select

    bd%alpha_grv = os_p%alpha_grv

    select case (os_p%time_factor)
    case ('OSC')
       bd%alpha_omg = 1._WP
    case ('EXP')
       bd%alpha_omg = -1._WP
    case default
       $ABORT(Invalid time_factor)
    end select

    call bd%stencil_(pt_i, pt_o)

    bd%n_i = 2
    bd%n_o = 2

    bd%n_e = 4

    ! Finish

    return
    
  end function ad_bound_t_

  !****

  subroutine stencil_ (this, pt_i, pt_o)

    class(ad_bound_t), intent(inout) :: this
    type(point_t), intent(in)        :: pt_i
    type(point_t), intent(in)        :: pt_o

    class(model_t), pointer :: ml

    ! Calculate coefficients at the stencil points

    ml => this%cx%model()

    call check_model(ml, [I_V_2,I_U,I_C_1,I_GAMMA_1])

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
    case (GAMMA_TYPE)
        this%coeff(2,J_C_1) = ml%coeff(I_C_1, pt_o)
        this%coeff(2,J_GAMMA_1) = ml%coeff(I_GAMMA_1, pt_o)
        this%coeff(2,J_V) = ml%coeff(I_V_2, pt_o)
    case default
       $ABORT(Invalid type_o)
    end select

    ! Set up stencil for the tr component

    call this%tr%stencil([pt_i,pt_o])

    ! Store the stencil points for on-the-fly evaluations

    this%pt = [pt_i,pt_o]

    ! Finish

    return

  end subroutine stencil_

  !****

  subroutine build_i (this, st, B, scl)

    class(ad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: B(:,:)
    real(WP), intent(out)         :: scl(:)

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

    class(ad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: B(:,:)
    real(WP), intent(out)         :: scl(:)

    real(WP) :: Omega_rot
    real(WP) :: omega_c
    real(WP) :: l_i

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
        
      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(l_i)
      B(2,4) = alpha_grv*(-1._WP) + (1._WP - alpha_grv)

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_regular_i_

  !****

  subroutine build_zero_r_i_ (this, st, B, scl)

    class(ad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: B(:,:)
    real(WP), intent(out)         :: scl(:)

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
        
      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(0._WP)
      B(2,4) = alpha_grv*(1._WP) + (1._WP - alpha_grv)

      scl = 1._WP

    end associate
      
    ! Finish

    return

  end subroutine build_zero_r_i_

  !****

  subroutine build_zero_h_i_ (this, st, B, scl)

    class(ad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: B(:,:)
    real(WP), intent(out)         :: scl(:)

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
        
      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(0._WP)
      B(2,4) = alpha_grv*(1._WP) + (1._WP - alpha_grv)

      scl = 1._WP

    end associate
      
    ! Finish

    return

  end subroutine build_zero_h_i_

  !****

  subroutine build_o (this, st, B, scl)

    class(ad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: B(:,:)
    real(WP), intent(out)         :: scl(:)

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
    case (GAMMA_TYPE)
       call this%build_gamma_o_(st, B, scl)
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

    class(ad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: B(:,:)
    real(WP), intent(out)         :: scl(:)

    real(WP) :: Omega_rot
    real(WP) :: l_e

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions (vacuum)

    associate( &
         U => this%coeff(2,J_U), &
         pt => this%pt(2), &
         alpha_grv => this%alpha_grv)

      Omega_rot = this%cx%Omega_rot(pt)

      l_e = this%cx%l_e(Omega_rot, st)

      ! Set up the boundary conditions

      B(1,1) = 1._WP
      B(1,2) = -1._WP
      B(1,3) = alpha_grv*(0._WP)
      B(1,4) = alpha_grv*(0._WP)
      
      B(2,1) = alpha_grv*(U)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(l_e + 1._WP) + (1._WP - alpha_grv)
      B(2,4) = alpha_grv*(1._WP)

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_vacuum_o_

  !****

  subroutine build_dziem_o_ (this, st, B, scl)

    class(ad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: B(:,:)
    real(WP), intent(out)         :: scl(:)

    real(WP) :: Omega_rot
    real(WP) :: omega_c
    real(WP) :: lambda
    real(WP) :: l_e

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions ([Dzi1971] formulation)

    associate( &
         V => this%coeff(2,J_V), &
         c_1 => this%coeff(2,J_C_1), &
         pt => this%pt(2), &
         alpha_grv => this%alpha_grv, &
         alpha_omg => this%alpha_omg)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, st)

      lambda = this%cx%lambda(Omega_rot, st)
      l_e = this%cx%l_e(Omega_rot, st)

      ! Set up the boundary conditions

      B(1,1) = 1._WP + (lambda/(c_1*alpha_omg*omega_c**2) - 4._WP - c_1*alpha_omg*omega_c**2)/V
      B(1,2) = -1._WP
      B(1,3) = alpha_grv*((lambda/(c_1*alpha_omg*omega_c**2) - l_e - 1._WP)/V)
      B(1,4) = alpha_grv*(0._WP)
      
      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(l_e + 1._WP) + (1._WP - alpha_grv)
      B(2,4) = alpha_grv*(1._WP)

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_dziem_o_

  !****

  subroutine build_decomp_o_ (this, st, B, scl)

    class(ad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: B(:,:)
    real(WP), intent(out)         :: scl(:)

    real(WP) :: Omega_rot
    real(WP) :: omega_c
    real(WP) :: lambda
    real(WP) :: l_e
    real(WP) :: a_11
    real(WP) :: a_12
    real(WP) :: a_13
    real(WP) :: a_14
    real(WP) :: a_21
    real(WP) :: a_22
    real(WP) :: a_23
    real(WP) :: a_24
    real(WP) :: chi
    real(WP) :: G_1
    real(WP) :: G_2

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
         pt => this%pt(2), &
         alpha_grv => this%alpha_grv, &
         alpha_omg => this%alpha_omg)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, st)

      lambda = this%cx%lambda(Omega_rot, st)
      l_e = this%cx%l_e(Omega_rot, st)
      
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

      chi = atmos_chi(V, As, c_1, Gamma_1, omega_c, lambda)

      ! Evaluate selected elements of the left eigenvectors

      G_1 = (l_e*H_(-l_e, chi) + (l_e+1._WP)*H_(l_e+1._WP, chi))/(2._WP*l_e + 1._WP)
      G_2 = (H_(-l_e, chi) - H_(l_e+1._WP, chi))/(2._WP*l_e + 1._WP)

      ! Set up the boundary conditions

      B(1,1) = -(chi - a_11)
      B(1,2) = a_12
      B(1,3) = -alpha_grv*G_1
      B(1,4) = alpha_grv*G_2

      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(l_e + 1._WP) + (1._WP - alpha_grv)
      B(2,4) = alpha_grv*(1._WP)

      scl = 1._WP

    end associate

    ! Finish

    return

  contains

    function H_ (l, chi) result (H)

      real(WP), intent(in) :: l
      real(WP), intent(in) :: chi
      real(WP)             :: H

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

    class(ad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: B(:,:)
    real(WP), intent(out)         :: scl(:)

    real(WP) :: Omega_rot
    real(WP) :: omega_c
    real(WP) :: lambda
    real(WP) :: l_e
    real(WP) :: chi
    real(WP) :: b_11
    real(WP) :: b_12

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions ([Chr2008] formulation)

    associate( &
         V => this%coeff(2,J_V), &
         As => this%coeff(2,J_AS), &
         c_1 => this%coeff(2,J_C_1), &
         Gamma_1 => this%coeff(2,J_GAMMA_1), &
         pt => this%pt(2), &
         alpha_grv => this%alpha_grv, &
         alpha_omg => this%alpha_omg)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, st)

      lambda = this%cx%lambda(Omega_rot, st)
      l_e = this%cx%l_e(Omega_rot, st)

      chi = atmos_chi(V, As, c_1, Gamma_1, omega_c, lambda)

      b_11 = V/Gamma_1 - 3._WP
      b_12 = lambda/(c_1*alpha_omg*omega_c**2) - V/Gamma_1

      ! Set up the boundary conditions

      B(1,1) = chi - b_11
      B(1,2) = -b_12
      B(1,3) = alpha_grv*((lambda/(c_1*alpha_omg*omega_c**2) - l_e - 1._WP)*b_12/(V/Gamma_1 + As))
      B(1,4) = alpha_grv*(0._WP)

      B(2,1) = alpha_grv*(0._WP)
      B(2,2) = alpha_grv*(0._WP)
      B(2,3) = alpha_grv*(l_e + 1._WP) + (1._WP - alpha_grv)
      B(2,4) = alpha_grv*(1._WP)

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_jcd_o_

  !****

  subroutine build_gamma_o_ (this, st, B, scl)

    class(ad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)  :: st
    real(WP), intent(out)         :: B(:,:)
    real(WP), intent(out)         :: scl(:)

    real(WP) :: Omega_rot
    real(WP) :: omega_c
    real(WP) :: l_e
    real(WP) :: lambda

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the inner boundary conditions (regular-enforcing)

    associate( &
         V => this%coeff(2,J_V), &
         c_1 => this%coeff(2,J_C_1), &
         Gamma_1 => this%coeff(2,J_GAMMA_1), &
         pt => this%pt(2), &
         alpha_grv => this%alpha_grv, &
         alpha_omg => this%alpha_omg)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, st)

      l_e = this%cx%l_e(Omega_rot, st)
      lambda = this%cx%lambda(Omega_rot, st)

      ! Set up the boundary conditions

      B(1,1) = 1._WP
      B(1,2) = 0._WP
      B(1,3) = 0._WP
      B(1,4) = 0._WP

      B(2,1) = 0._WP
      B(2,2) = 1._WP !lambda/(c_1*alpha_omg*omega_c**2) !- V/Gamma_1 * alpha_gamma
      B(2,3) = alpha_grv * 1._WP !alpha_grv*lambda/(c_1*alpha_omg*omega_c**2)
      B(2,4) = 0._WP

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_gamma_o_

end module gyre_ad_bound
