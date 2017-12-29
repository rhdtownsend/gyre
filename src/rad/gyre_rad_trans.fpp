! Module   : gyre_rad_vars
! Purpose  : adiabatic radial variables transformations
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

module gyre_rad_trans

  ! Uses

  use core_kinds

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

  integer, parameter :: GYRE_SET = 0
  integer, parameter :: DZIEM_SET = 1
  integer, parameter :: JCD_SET = 2
  integer, parameter :: MIX_SET = 3
  integer, parameter :: LAGP_SET = 4

  integer, parameter :: J_V_2 = 1
  integer, parameter :: J_DV_2 = 2
  integer, parameter :: J_C_1 = 3
  integer, parameter :: J_DC_1 = 4

  integer, parameter :: J_LAST = J_DC_1

  ! Derived-type definitions

  type :: rad_trans_t
     private
     type(context_t), pointer :: cx => null()
     real(WP), allocatable    :: coeff(:,:)
     integer                  :: set
     integer                  :: n_e
   contains
     private
     procedure, public :: stencil
     procedure, public :: trans_eqns
     procedure, public :: trans_cond
     procedure, public :: trans_vars
     procedure         :: G_
     procedure         :: G_jcd_
     procedure         :: G_lagp_
     procedure         :: H_
     procedure         :: H_jcd_
     procedure         :: H_lagp_
     procedure         :: dH_
     procedure         :: dH_jcd_
     procedure         :: dH_lagp_
  end type rad_trans_t

  ! Interfaces

  interface rad_trans_t
     module procedure rad_trans_t_
  end interface rad_trans_t

  ! Access specifiers

  private

  public :: rad_trans_t

  ! Procedures

contains

  function rad_trans_t_ (cx, pt_i, md_p, os_p) result (tr)

    type(context_t), pointer, intent(in) :: cx
    type(point_t), intent(in)            :: pt_i
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(rad_trans_t)                    :: tr

    ! Construct the rad_trans_t

    tr%cx => cx

    select case (os_p%variables_set)
    case ('GYRE')
       tr%set = GYRE_SET
    case ('DZIEM')
       tr%set = DZIEM_SET
    case ('JCD')
       tr%set = JCD_SET
    case ('MIX')
       tr%set = MIX_SET
    case ('LAGP')
       tr%set = LAGP_SET
    case default
       $ABORT(Invalid variables_set)
    end select

    tr%n_e = 2

    ! Finish

    return

  end function rad_trans_t_

  !****

  subroutine stencil (this, pt)

    class(rad_trans_t), intent(inout) :: this
    type(point_t), intent(in)         :: pt(:)

    integer :: n_s
    integer :: i

    ! Calculate coefficients at the stencil points

    associate (ml => this%cx%ml)

      call check_model(ml, [I_V_2,I_U,I_C_1])

      n_s = SIZE(pt)

      if (ALLOCATED(this%coeff)) deallocate(this%coeff)
      allocate(this%coeff(n_s,J_LAST))

      do i = 1, n_s
         if (ml%is_vacuum(pt(i))) then
            this%coeff(i,J_V_2) = HUGE(0._WP)
            this%coeff(i,J_DV_2) = HUGE(0._WP)
         else
            this%coeff(i,J_V_2) = ml%coeff(I_V_2, pt(i))
            this%coeff(i,J_DV_2) = ml%dcoeff(I_V_2, pt(i))
         endif
         this%coeff(i,J_C_1) = ml%coeff(I_C_1, pt(i))
         this%coeff(i,J_DC_1) = ml%dcoeff(I_C_1, pt(i))
      end do

    end associate

    ! Finish

    return

  end subroutine stencil

  !****

  subroutine trans_eqns (this, xA, i, st, from)

    class(rad_trans_t), intent(in) :: this
    real(WP), intent(inout)        :: xA(:,:)
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    logical, intent(in), optional  :: from

    logical  :: from_
    real(WP) :: G(this%n_e,this%n_e)
    real(WP) :: H(this%n_e,this%n_e)
    real(WP) :: dH(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(xA, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(xA, 2),this%n_e)

    if (PRESENT(from)) then
       from_ = from
    else
       from_ = .TRUE.
    endif

    ! Transform equations to/from GYRE's canonical form

    if (from_) then

       ! Convert from

       if (this%set /= GYRE_SET .AND. &
           this%set /= DZIEM_SET .AND. &
           this%set /= MIX_SET) then
          G = this%G_(i, st)
          H = this%H_(i, st)
          dH = this%dH_(i, st)
          xA = MATMUL(G, MATMUL(xA, H) - dH)
       endif

    else

       ! Convert to

       $ABORT(Not currently supported)

    end if

    ! Finish

    return

  end subroutine trans_eqns
  
  !****

  subroutine trans_cond (this, C, i, st, from)

    class(rad_trans_t), intent(in) :: this
    real(WP), intent(inout)        :: C(:,:)
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    logical, intent(in), optional  :: from

    logical  :: from_
    real(WP) :: G(this%n_e,this%n_e)
    real(WP) :: H(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(C, 2),this%n_e)

    if (PRESENT(from)) then
       from_ = from
    else
       from_ = .TRUE.
    endif

    ! Transform boundary/match conditions to/from GYRE's canonical form

    if (from_) then

       ! Convert from

       if (this%set /= GYRE_SET .AND. &
           this%set /= DZIEM_SET .AND. &
           this%set /= MIX_SET) then
          H = this%H_(i, st)
          C = MATMUL(C, H)
       endif

    else

       ! Convert to

       if (this%set /= GYRE_SET .AND. &
           this%set /= DZIEM_SET .AND. &
           this%set /= MIX_SET) then
          G = this%G_(i, st)
          C = MATMUL(C, G)
       endif

    end if

    ! Finish

    return

  end subroutine trans_cond

  !****

  subroutine trans_vars (this, y, i, st, from)

    class(rad_trans_t), intent(in) :: this
    real(WP), intent(inout)        :: y(:)
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    logical, intent(in), optional  :: from

    logical  :: from_
    real(WP) :: G(this%n_e,this%n_e)
    real(WP) :: H(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(y),this%n_e)

    if (PRESENT(from)) then
       from_ = from
    else
       from_ = .TRUE.
    endif

    ! Convert variables to/from GYRE's canonical form

    if (from_) then

       ! Convert from

       if (this%set /= GYRE_SET .AND. &
           this%set /= DZIEM_SET .AND. &
           this%set /= MIX_SET) then
          G = this%G_(i, st)
          y = MATMUL(G, y)
       endif

    else

       ! Convert to

       if (this%set /= GYRE_SET .AND. &
           this%set /= DZIEM_SET .AND. &
           this%set /= MIX_SET) then
          H = this%H_(i, st)
          y = MATMUL(H, y)
       endif

    end if

    ! Finish

    return

  end subroutine trans_vars

  !****

  function G_ (this, i, st) result (G)

    class(rad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: G(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert variables from
    ! GYRE's canonical form

    select case (this%set)
    case (JCD_SET)
       G = this%G_jcd_(i, st)
    case (LAGP_SET)
       G = this%G_lagp_(i, st)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return
    
  end function G_
  
  !****

  function G_jcd_ (this, i, st) result (G)

    class(rad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: G(this%n_e,this%n_e)

    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! from GYRE's canonical form

    associate( &
         omega => st%omega, &
         c_1 => this%coeff(i,J_C_1))

      omega_c = omega

      ! Set up the matrix

      G(1,1) = 1._WP
      G(1,2) = 0._WP

      G(2,1) = 0._WP
      G(2,2) = 1._WP/(c_1*omega_c**2)

    end associate
      
    ! Finish

    return

  end function G_jcd_

  !****

  function G_lagp_ (this, i, st) result (G)

    class(rad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: G(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert LAGP variables
    ! from GYRE's canonical form

    associate( &
         V_2 => this%coeff(i,J_V_2))

      ! Set up the matrix

      G(1,1) = 1._WP
      G(1,2) = 0._WP

      G(2,1) = -V_2
      G(2,2) = V_2

    end associate

    ! Finish

    return

  end function G_lagp_

  !****

  function H_ (this, i, st) result (H)

    class(rad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: H(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert variables to
    ! canonical form

    select case (this%set)
    case (JCD_SET)
       H = this%H_jcd_(i, st)
    case (LAGP_SET)
       H = this%H_lagp_(i, st)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function H_

  !****

  function H_jcd_ (this, i, st) result (H)

    class(rad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: H(this%n_e,this%n_e)

    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! to GYRE's canonical form

    associate( &
         omega => st%omega, &
         c_1 => this%coeff(i,J_C_1))

      omega_c = omega

      ! Set up the matrix
      
      H(1,1) = 1._WP
      H(1,2) = 0._WP
      
      H(2,1) = 0._WP
      H(2,2) = c_1*omega_c**2

    end associate

    ! Finish

    return

  end function H_jcd_

  !****

  function H_lagp_ (this, i, st) result (H)

    class(rad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: H(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert LAGP variables
    ! to GYRE's canonical form

    associate( &
         V_2 => this%coeff(i,J_V_2))

      ! Set up the matrix

      H(1,1) = 1._WP
      H(1,2) = 0._WP
      
      H(2,1) = 1._WP
      H(2,2) = 1._WP/V_2

    end associate

    ! Finish

    return

  end function H_lagp_

  !****

  function dH_ (this, i, st) result (dH)

    class(rad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: dH(this%n_e,this%n_e)

    ! Evaluate the derivative x dH/dx of the transformation matrix H

    select case (this%set)
    case (JCD_SET)
       dH = this%dH_jcd_(i, st)
    case (LAGP_SET)
       dH = this%dH_lagp_(i, st)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function dH_

  !****

  function dH_jcd_ (this, i, st) result (dH)

    class(rad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: dH(this%n_e,this%n_e)

    real(WP) :: omega_c

    ! Evaluate the derivative x dH/dx of the JCD-variables
    ! transformation matrix H

    associate( &
         omega => st%omega, &
         c_1 => this%coeff(i,J_C_1), &
         dc_1 => this%coeff(i,J_DC_1))

      omega_c = omega

      ! Set up the matrix
      
      dH(1,1) = 0._WP
      dH(1,2) = 0._WP

      dH(2,1) = 0._WP
      dH(2,2) = c_1*dc_1*omega_c**2

    end associate

    ! Finish

    return

  end function dH_jcd_

  !****

  function dH_lagp_ (this, i, st) result (dH)

    class(rad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(r_state_t), intent(in)   :: st
    real(WP)                       :: dH(this%n_e,this%n_e)

    ! Evaluate the derivative x dH/dx of the LAGP-variables
    ! transformation matrix H

    associate( &
         V_2 => this%coeff(i,J_V_2), &
         dV_2 => this%coeff(i,J_DV_2))

      ! Set up the matrix

      dH(1,1) = 0._WP
      dH(1,2) = 0._WP

      dH(2,1) = 0._WP
      dH(2,2) = -dV_2/V_2

    end associate

    ! Finish

    return

  end function dH_lagp_
  
end module gyre_rad_trans
