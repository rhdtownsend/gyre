! Incfile  : gyre_rad_bound
! Purpose  : adiabatic radial boundary conditions
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

module gyre_rad_bound

  ! Uses

  use core_kinds

  use gyre_atmos
  use gyre_bound
  use gyre_context
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rad_trans
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: REGULAR_TYPE = 1
  integer, parameter :: ZERO_R_TYPE = 2
  integer, parameter :: VACUUM_TYPE = 3
  integer, parameter :: DZIEM_TYPE = 4
  integer, parameter :: UNNO_TYPE = 5
  integer, parameter :: JCD_TYPE = 6
  integer, parameter :: LUAN_TYPE = 7

  integer, parameter :: J_V = 1
  integer, parameter :: J_V_G = 2
  integer, parameter :: J_AS = 3
  integer, parameter :: J_U = 4
  integer, parameter :: J_C_1 = 5

  integer, parameter :: J_LAST = J_C_1

  ! Derived-type definitions

  type, extends (r_bound_t) :: rad_bound_t
     private
     type(context_t), pointer :: cx => null()
     type(rad_trans_t)        :: tr
     real(WP), allocatable    :: coeff(:,:)
     real(WP)                 :: alpha_om
     integer                  :: type_i
     integer                  :: type_o
   contains 
     private
     procedure         :: stencil_
     procedure, public :: build_i
     procedure         :: build_regular_i_
     procedure         :: build_zero_r_i_
     procedure, public :: build_o
     procedure         :: build_vacuum_o_
     procedure         :: build_dziem_o_
     procedure         :: build_unno_o_
     procedure         :: build_jcd_o_
     procedure         :: build_luan_o_
  end type rad_bound_t

  ! Interfaces

  interface rad_bound_t
     module procedure rad_bound_t_
  end interface rad_bound_t

  ! Access specifiers

  private

  public :: rad_bound_t

  ! Procedures

contains

  function rad_bound_t_ (cx, pt_i, pt_o, md_p, os_p) result (bd)

    type(context_t), pointer, intent(in) :: cx
    type(point_t), intent(in)            :: pt_i
    type(point_t), intent(in)            :: pt_o
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(rad_bound_t)                    :: bd

    ! Construct the ad_bound_t

    bd%cx => cx
    
    bd%tr = rad_trans_t(cx, pt_i, md_p, os_p)

    select case (os_p%inner_bound)
    case ('REGULAR')
       $ASSERT(pt_i%x == 0._WP,Boundary condition invalid for x /= 0)
       bd%type_i = REGULAR_TYPE
    case ('ZERO_R')
       $ASSERT(pt_i%x /= 0._WP,Boundary condition invalid for x == 0)
       bd%type_i = ZERO_R_TYPE
    case default
       $ABORT(Invalid inner_bound)
    end select

    select case (os_p%outer_bound)
    case ('VACUUM')
       bd%type_o = VACUUM_TYPE
    case ('DZIEM')
       if (cx%ml%is_vacuum(pt_o)) then
          bd%type_o = VACUUM_TYPE
       else
          bd%type_o = DZIEM_TYPE
       endif
    case ('UNNO')
       if (cx%ml%is_vacuum(pt_o)) then
          bd%type_o = VACUUM_TYPE
       else
          bd%type_o = UNNO_TYPE
       endif
    case ('JCD')
       if (cx%ml%is_vacuum(pt_o)) then
          bd%type_o = VACUUM_TYPE
       else
          bd%type_o = JCD_TYPE
       endif
    case ('LUAN')
       if (cx%ml%is_vacuum(pt_o)) then
          bd%type_o = VACUUM_TYPE
       else
          bd%type_o = LUAN_TYPE
       endif
    case default
       $ABORT(Invalid outer_bound)
    end select

    select case (os_p%time_factor)
    case ('OSC')
       bd%alpha_om = 1._WP
    case ('EXP')
       bd%alpha_om = -1._WP
    case default
       $ABORT(Invalid time_factor)
    end select

    call bd%stencil_(pt_i, pt_o)

    bd%n_i = 1
    bd%n_o = 1

    bd%n_e = 2

    ! Finish

    return
    
  end function rad_bound_t_

  !****

  subroutine stencil_ (this, pt_i, pt_o)

    class(rad_bound_t), intent(inout) :: this
    type(point_t), intent(in)         :: pt_i
    type(point_t), intent(in)         :: pt_o

    ! Calculate coefficients at the stencil points

    associate (ml => this%cx%ml)

      call check_model(ml, [I_V_2,I_U,I_C_1])

      allocate(this%coeff(2,J_LAST))

      ! Inner boundary

      select case (this%type_i)
      case (REGULAR_TYPE)
         this%coeff(1,J_C_1) = ml%coeff(I_C_1, pt_i)
      case (ZERO_R_TYPE)
      case default
         $ABORT(Invalid type_i)
      end select

      ! Outer boundary

      select case (this%type_o)
      case (VACUUM_TYPE)
         this%coeff(2, J_U) = ml%coeff(I_U, pt_o)
      case (DZIEM_TYPE)
         this%coeff(2,J_V) = ml%coeff(I_V_2, pt_o)*pt_o%x**2
         this%coeff(2,J_C_1) = ml%coeff(I_C_1, pt_o)
      case (UNNO_TYPE)
         call eval_atmos_coeffs_unno(ml, pt_o, this%coeff(2,J_V_G), &
              this%coeff(2,J_AS), this%coeff(2,J_U), this%coeff(2,J_C_1))
      case (JCD_TYPE)
         call eval_atmos_coeffs_jcd(ml, pt_o, this%coeff(2,J_V_G), &
              this%coeff(2,J_AS), this%coeff(2,J_U), this%coeff(2,J_C_1))
      case (LUAN_TYPE)
         call eval_atmos_coeffs_luan(ml, pt_o, this%coeff(2,J_V_G), &
              this%coeff(2,J_AS), this%coeff(2,J_U), this%coeff(2,J_C_1))
      case default
         $ABORT(Invalid type_o)
      end select

    end associate

    ! Set up stencil for the tr component

    call this%tr%stencil([pt_i,pt_o])

    ! Finish

    return

  end subroutine stencil_

  !****

  subroutine build_i (this, st, B, scl)

    class(rad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_i)
    
    ! Evaluate the inner boundary conditions

    select case (this%type_i)
    case (REGULAR_TYPE)
       call this%build_regular_i_(st, B, scl)
    case (ZERO_R_TYPE)
       call this%build_zero_r_i_(st, B, scl)
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

    class(rad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    real(WP) :: omega_c

    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(scl),this%n_i)

    ! Evaluate the inner boundary conditions (regular-enforcing)

    associate( &
         omega => st%omega, &
         c_1 => this%coeff(1,J_C_1), &
         alpha_om => this%alpha_om)

      omega_c = omega

      ! Set up the boundary conditions

      B(1,1) = c_1*alpha_om*omega_c**2
      B(1,2) = 0._WP

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_regular_i_

  !****

  subroutine build_zero_r_i_ (this, st, B, scl)

    class(rad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_i)

    ! Evaluate the inner boundary conditions (zero radial
    ! displacement)

    ! Set up the boundary conditions

    B(1,1) = 1._WP
    B(1,2) = 0._WP

    scl = 1._WP
      
    ! Finish

    return

  end subroutine build_zero_r_i_

  !****

  subroutine build_o (this, st, B, scl)

    class(rad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions

    select case (this%type_o)
    case (VACUUM_TYPE)
       call this%build_vacuum_o_(st, B, scl)
    case (DZIEM_TYPE)
       call this%build_dziem_o_(st, B, scl)
    case (UNNO_TYPE)
       call this%build_unno_o_(st, B, scl)
    case (JCD_TYPE)
       call this%build_jcd_o_(st, B, scl)
    case (LUAN_TYPE)
       call this%build_luan_o_(st, B, scl)
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

    class(rad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions (vacuum)

    ! Set up the boundary conditions

    B(1,1) = 1._WP
    B(1,2) = -1._WP
      
    scl = 1._WP
    
    ! Finish

    return

  end subroutine build_vacuum_o_

  !****

  subroutine build_dziem_o_ (this, st, B, scl)

    class(rad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    real(WP) :: omega_c

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions ([Dzi1971] formulation)

    associate( &
         omega => st%omega, &
         V => this%coeff(2,J_V), &
         c_1 => this%coeff(2,J_C_1), &
         alpha_om => this%alpha_om)

      omega_c = omega

      ! Set up the boundary conditions
        
      B(1,1) = 1 - (4._WP + c_1*alpha_om*omega_c**2)/V
      B(1,2) = -1._WP
      
      scl = 1._WP
         
    end associate

    ! Finish

    return

  end subroutine build_dziem_o_

  !****
  
  subroutine build_unno_o_ (this, st, B, scl)

    class(rad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    real(WP) :: omega_c
    real(WP) :: beta
    real(WP) :: b_11
    real(WP) :: b_12

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions ([Unn1989] formulation)

    associate( &
         omega => st%omega, &
         V_g => this%coeff(2,J_V_G), &
         As => this%coeff(2,J_AS), &
         U => this%coeff(2,J_U), &
         c_1 => this%coeff(2,J_C_1))

      omega_c = omega

      beta = atmos_beta(V_g, As, U, c_1, omega_c, 0._WP)
      
      b_11 = V_g - 3._WP
      b_12 = -V_g
    
      ! Set up the boundary conditions
      
      B(1,1) = beta - b_11
      B(1,2) = -b_12

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_unno_o_

  !****

  subroutine build_jcd_o_ (this, st, B, scl)

    class(rad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    real(WP) :: omega_c
    real(WP) :: beta
    real(WP) :: b_11
    real(WP) :: b_12

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions ([Chr2008] formulation)

    ! Calculate coefficients

    associate( &
         omega => st%omega, &
         V_g => this%coeff(2,J_V_G), &
         As => this%coeff(2,J_AS), &
         U => this%coeff(2,J_U), &
         c_1 => this%coeff(2,J_C_1))

      omega_c = omega

      beta = atmos_beta(V_g, As, U, c_1, omega_c, 0._WP)

      b_11 = V_g - 3._WP
      b_12 = -V_g

      ! Set up the boundary conditions

      B(1,1) = beta - b_11
      B(1,2) = -b_12

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_jcd_o_

  !****

  subroutine build_luan_o_ (this, st, B, scl)

    class(rad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    real(WP) :: omega_c
    real(WP) :: beta
    real(WP) :: b_11
    real(WP) :: b_12

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions (Luan formulation)

    ! Calculate coefficients

    associate( &
         omega => st%omega, &
         V_g => this%coeff(2,J_V_G), &
         As => this%coeff(2,J_AS), &
         U => this%coeff(2,J_U), &
         c_1 => this%coeff(2,J_C_1))

      omega_c = omega

      beta = atmos_beta(V_g, As, U, c_1, omega_c, 0._WP)

      b_11 = V_g - 3._WP
      b_12 = -V_g

      ! Set up the boundary conditions

      B(1,1) = beta - b_11
      B(1,2) = -b_12

      scl = 1._WP

    end associate

    ! Finish

    return

  end subroutine build_luan_o_

end module gyre_rad_bound
