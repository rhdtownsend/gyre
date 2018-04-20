! Module   : gyre_ad_vars
! Purpose  : nonadiabatic variables transformations
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

module gyre_nad_trans

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
  integer, parameter :: LAGP_SET = 3

  integer, parameter :: J_V_2 = 1
  integer, parameter :: J_DV_2 = 2
  integer, parameter :: J_U = 3
  integer, parameter :: J_DU = 4
  integer, parameter :: J_C_1 = 5
  integer, parameter :: J_DC_1 = 6
  integer, parameter :: J_OMEGA_ROT = 7

  integer, parameter :: J_LAST = J_OMEGA_ROT

  ! Derived-type definitions

  type :: nad_trans_t
     private
     type(context_t), pointer :: cx => null()
     real(WP), allocatable    :: coeff(:,:)
     integer                  :: set
     integer                  :: l
     integer                  :: m
     integer                  :: n_e
   contains
     private
     procedure, public :: stencil
     procedure, public :: trans_eqns
     procedure, public :: trans_cond
     procedure, public :: trans_vars
     procedure         :: G_
     procedure         :: G_dziem_
     procedure         :: G_jcd_
     procedure         :: G_lagp_
     procedure         :: H_
     procedure         :: H_dziem_
     procedure         :: H_jcd_
     procedure         :: H_lagp_
     procedure         :: dH_
     procedure         :: dH_jcd_
     procedure         :: dH_lagp_
  end type nad_trans_t
  
  ! Interfaces

  interface nad_trans_t
     module procedure nad_trans_t_
  end interface nad_trans_t

  ! Access specifiers

  private

  public :: nad_trans_t

  ! Procedures

contains

  function nad_trans_t_ (cx, md_p, os_p) result (tr)

    type(context_t), pointer, intent(in) :: cx
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(nad_trans_t)                    :: tr

    ! Construct the nad_trans_t

    tr%cx => cx

    select case (os_p%variables_set)
    case ('GYRE')
       tr%set = GYRE_SET
    case ('DZIEM')
       tr%set = DZIEM_SET
    case ('JCD')
       tr%set = JCD_SET
    case ('LAGP')
       tr%set = LAGP_SET
    case default
       $ABORT(Invalid variables_set)
    end select

    tr%l = md_p%l
    tr%m = md_p%m

    tr%n_e = 6

    ! Finish

    return

  end function nad_trans_t_

  !****

  subroutine stencil (this, pt)

    class(nad_trans_t), intent(inout) :: this
    type(point_t), intent(in)         :: pt(:)

    integer :: n_s
    integer :: i

    ! Calculate coefficients at the stencil points

    associate (ml => this%cx%ml)

      call check_model(ml, [I_V_2,I_U,I_C_1,I_OMEGA_ROT])

      n_s = SIZE(pt)

      if (ALLOCATED(this%coeff)) deallocate(this%coeff)
      allocate(this%coeff(n_s,J_LAST))

      do i = 1, n_s
         if (ml%is_vacuum(pt(i))) then
            this%coeff(i,J_V_2) = HUGE(0._WP)
            this%coeff(i,J_DV_2) = HUGE(0._WP)
            this%coeff(i,J_DU) = -HUGE(0._WP)
         else
            this%coeff(i,J_V_2) = ml%coeff(I_V_2, pt(i))
            this%coeff(i,J_DV_2) = ml%dcoeff(I_V_2, pt(i))
            this%coeff(i,J_DU) = ml%dcoeff(I_U, pt(i))
         endif
         this%coeff(i,J_U) = ml%coeff(I_U, pt(i))
         this%coeff(i,J_C_1) = ml%coeff(I_C_1, pt(i))
         this%coeff(i,J_DC_1) = ml%dcoeff(I_C_1, pt(i))
         this%coeff(i,J_OMEGA_ROT) = ml%coeff(I_OMEGA_ROT, pt(i))
      end do

    end associate

    ! Finish

    return

  end subroutine stencil

  !****

  subroutine trans_eqns (this, xA, i, st, from)

    class(nad_trans_t), intent(in) :: this
    complex(WP), intent(inout)     :: xA(:,:)
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    logical, intent(in), optional  :: from

    logical     :: from_
    complex(WP) :: G(this%n_e,this%n_e)
    complex(WP) :: H(this%n_e,this%n_e)
    complex(WP) :: dH(this%n_e,this%n_e)

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

       if (this%set /= GYRE_SET) then
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

    class(nad_trans_t), intent(in) :: this
    complex(WP), intent(inout)     :: C(:,:)
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    logical, intent(in), optional  :: from

    logical     :: from_
    complex(WP) :: G(this%n_e,this%n_e)
    complex(WP) :: H(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(C, 2),this%n_e)

    if (PRESENT(from)) then
       from_ = from
    else
       from_ = .TRUE.
    endif

    ! Transform boundary/match conditions to/from GYRE's canonical form

    if (from_) then

       ! Convert from

       if (this%set /= GYRE_SET) then
          H = this%H_(i, st)
          C = MATMUL(C, H)
       endif

    else

       ! Convert to

       if (this%set /= GYRE_SET) then
          G = this%G_(i, st)
          C = MATMUL(C, G)
       endif

    end if

    ! Finish

    return

  end subroutine trans_cond

  !****

  subroutine trans_vars (this, y, i, st, from)

    class(nad_trans_t), intent(in) :: this
    complex(WP), intent(inout)     :: y(:)
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    logical, intent(in), optional  :: from

    logical     :: from_
    complex(WP) :: G(this%n_e,this%n_e)
    complex(WP) :: H(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(y),this%n_e)

    if (PRESENT(from)) then
       from_ = from
    else
       from_ = .TRUE.
    endif

    ! Convert variables to/from GYRE's canonical form

    if (from_) then

       ! Convert from

       if (this%set /= GYRE_SET) then
          G = this%G_(i, st)
          y = MATMUL(G, y)
       endif

    else

       ! Convert to

       if (this%set /= GYRE_SET) then
          H = this%H_(i, st)
          y = MATMUL(H, y)
       endif

    end if

    ! Finish

    return

  end subroutine trans_vars

  !****

  function G_ (this, i, st) result (G)

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: G(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert variables from
    ! the canonical form

    select case (this%set)
    case (DZIEM_SET)
       G = this%G_dziem_(i, st)
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

  function G_dziem_ (this, i, st) result (G)

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: G(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert DZIEM variables
    ! from the canonical form

    ! Set up the matrix

    G(1,1) = 1._WP
    G(1,2) = 0._WP
    G(1,3) = 0._WP
    G(1,4) = 0._WP
    G(1,5) = 0._WP
    G(1,6) = 0._WP

    G(2,1) = 0._WP
    G(2,2) = 1._WP
    G(2,3) = 1._WP
    G(2,4) = 0._WP
    G(2,5) = 0._WP
    G(2,6) = 0._WP

    G(3,1) = 0._WP
    G(3,2) = 0._WP
    G(3,3) = 1._WP
    G(3,4) = 0._WP
    G(3,5) = 0._WP
    G(3,6) = 0._WP

    G(4,1) = 0._WP
    G(4,2) = 0._WP
    G(4,3) = 0._WP
    G(4,4) = 1._WP
    G(4,5) = 0._WP
    G(4,6) = 0._WP

    G(5,1) = 0._WP
    G(5,2) = 0._WP
    G(5,3) = 0._WP
    G(5,4) = 0._WP
    G(5,5) = 1._WP
    G(5,6) = 0._WP

    G(6,1) = 0._WP
    G(6,2) = 0._WP
    G(6,3) = 0._WP
    G(6,4) = 0._WP
    G(6,5) = 0._WP
    G(6,6) = 1._WP

    ! Finish

    return

  end function G_dziem_

  !****

  function G_jcd_ (this, i, st) result (G)

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: G(this%n_e,this%n_e)

    complex(WP) :: omega_c
    complex(WP) :: lambda

    ! Evaluate the transformation matrix to convert JCD variables
    ! from the canonical form

    associate ( &
         U => this%coeff(i,J_U), &
         c_1 => this%coeff(i,J_C_1), &
         Omega_rot => this%coeff(i,J_OMEGA_ROT))

      omega_c = this%cx%omega_c(Omega_rot, st)

      lambda = this%cx%lambda(Omega_rot, st)

      ! Set up the matrix

      if (this%l /= 0) then

         G(1,1) = 1._WP
         G(1,2) = 0._WP
         G(1,3) = 0._WP
         G(1,4) = 0._WP
         G(1,5) = 0._WP
         G(1,6) = 0._WP

         G(2,1) = 0._WP
         G(2,2) = lambda/(c_1*omega_c**2)
         G(2,3) = lambda/(c_1*omega_c**2)
         G(2,4) = 0._WP
         G(2,5) = 0._WP
         G(2,6) = 0._WP
         
         G(3,1) = 0._WP
         G(3,2) = 0._WP
         G(3,3) = -1._WP
         G(3,4) = 0._WP
         G(3,5) = 0._WP
         G(3,6) = 0._WP

         G(4,1) = 0._WP
         G(4,2) = 0._WP
         G(4,3) = -(1._WP - U)
         G(4,4) = -1._WP
         G(4,5) = 0._WP
         G(4,6) = 0._WP

         G(5,1) = 0._WP
         G(5,2) = 0._WP
         G(5,3) = 0._WP
         G(5,4) = 0._WP
         G(5,5) = 1._WP
         G(5,6) = 0._WP

         G(6,1) = 0._WP
         G(6,2) = 0._WP
         G(6,3) = 0._WP
         G(6,4) = 0._WP
         G(6,5) = 0._WP
         G(6,6) = 1._WP
          
      else

         G(1,1) = 1._WP
         G(1,2) = 0._WP
         G(1,3) = 0._WP
         G(1,4) = 0._WP
         G(1,5) = 0._WP
         G(1,6) = 0._WP

         G(2,1) = 0._WP
         G(2,2) = 1._WP/(c_1*omega_c**2)
         G(2,3) = 1._WP/(c_1*omega_c**2)
         G(2,4) = 0._WP
         G(2,5) = 0._WP
         G(2,6) = 0._WP

         G(3,1) = 0._WP
         G(3,2) = 0._WP
         G(3,3) = -1._WP
         G(3,4) = 0._WP
         G(3,5) = 0._WP
         G(3,6) = 0._WP

         G(4,1) = 0._WP
         G(4,2) = 0._WP
         G(4,3) = -(1._WP - U)
         G(4,4) = -1._WP
         G(4,5) = 0._WP
         G(4,6) = 0._WP

         G(5,1) = 0._WP
         G(5,2) = 0._WP
         G(5,3) = 0._WP
         G(5,4) = 0._WP
         G(5,5) = 1._WP
         G(5,6) = 0._WP

         G(6,1) = 0._WP
         G(6,2) = 0._WP
         G(6,3) = 0._WP
         G(6,4) = 0._WP
         G(6,5) = 0._WP
         G(6,6) = 1._WP

      end if

    end associate

    ! Finish

    return

  end function G_jcd_

  !****

  function G_lagp_ (this, i, st) result (G)

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: G(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert LAGP variables
    ! from the canonical form

    associate ( &
         V_2 => this%coeff(i,J_V_2))

      ! Set up the matrix
      
      G(1,1) = 1._WP
      G(1,2) = 0._WP
      G(1,3) = 0._WP
      G(1,4) = 0._WP
      G(1,5) = 0._WP
      G(1,6) = 0._WP

      G(2,1) = -V_2
      G(2,2) = V_2
      G(2,3) = 0._WP
      G(2,4) = 0._WP
      G(2,5) = 0._WP
      G(2,6) = 0._WP

      G(3,1) = 0._WP
      G(3,2) = 0._WP
      G(3,3) = 1._WP
      G(3,4) = 0._WP
      G(3,5) = 0._WP
      G(3,6) = 0._WP

      G(4,1) = 0._WP
      G(4,2) = 0._WP
      G(4,3) = 0._WP
      G(4,4) = 1._WP
      G(4,5) = 0._WP
      G(4,6) = 0._WP

      G(5,1) = 0._WP
      G(5,2) = 0._WP
      G(5,3) = 0._WP
      G(5,4) = 0._WP
      G(5,5) = 1._WP
      G(5,6) = 0._WP

      G(6,1) = 0._WP
      G(6,2) = 0._WP
      G(6,3) = 0._WP
      G(6,4) = 0._WP
      G(6,5) = 0._WP
      G(6,6) = 1._WP

    end associate

    ! Finish

    return

  end function G_lagp_

  !****

  function H_ (this, i, st) result (H)

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: H(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert variables to
    ! canonical form

    select case (this%set)
    case (DZIEM_SET)
       H = this%H_dziem_(i, st)
    case (JCD_SET)
       H = this%H_jcd_(i, st)
    case (LAGP_SET)
       H = this%H_lagp_(i, st)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function H_

  !****

  function H_dziem_ (this, i, st) result (H)

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: H(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert DZIEM variables
    ! to the canonical form

    ! Set up the matrix

    H(1,1) = 1._WP
    H(1,2) = 0._WP
    H(1,3) = 0._WP
    H(1,4) = 0._WP
    H(1,5) = 0._WP
    H(1,6) = 0._WP

    H(2,1) = 0._WP
    H(2,2) = 1._WP
    H(2,3) = -1._WP
    H(2,4) = 0._WP
    H(2,5) = 0._WP
    H(2,6) = 0._WP

    H(3,1) = 0._WP
    H(3,2) = 0._WP
    H(3,3) = 1._WP
    H(3,4) = 0._WP
    H(3,5) = 0._WP
    H(3,6) = 0._WP

    H(4,1) = 0._WP
    H(4,2) = 0._WP
    H(4,3) = 0._WP
    H(4,4) = 1._WP
    H(4,5) = 0._WP
    H(4,6) = 0._WP
         
    H(5,1) = 0._WP
    H(5,2) = 0._WP
    H(5,3) = 0._WP
    H(5,4) = 0._WP
    H(5,5) = 1._WP
    H(5,6) = 0._WP

    H(6,1) = 0._WP
    H(6,2) = 0._WP
    H(6,3) = 0._WP
    H(6,4) = 0._WP
    H(6,5) = 0._WP
    H(6,6) = 1._WP
    
    ! Finish

    return

  end function H_dziem_

  !****

  function H_jcd_ (this, i, st) result (H)

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: H(this%n_e,this%n_e)

    complex(WP) :: omega_c
    complex(WP) :: lambda

    ! Evaluate the transformation matrix to convert JCD variables
    ! to the canonical form

    associate ( &
         U => this%coeff(i,J_U), &
         c_1 => this%coeff(i,J_C_1), &
         Omega_rot => this%coeff(i,J_OMEGA_ROT))

      omega_c = this%cx%omega_c(Omega_rot, st)

      lambda = this%cx%lambda(Omega_rot, st)

      ! Set up the matrix
      
      if (this%l /= 0) then

         H(1,1) = 1._WP
         H(1,2) = 0._WP
         H(1,3) = 0._WP
         H(1,4) = 0._WP
         H(1,5) = 0._WP
         H(1,6) = 0._WP

         H(2,1) = 0._WP
         H(2,2) = c_1*omega_c**2/lambda
         H(2,3) = 1._WP
         H(2,4) = 0._WP
         H(2,5) = 0._WP
         H(2,6) = 0._WP

         H(3,1) = 0._WP
         H(3,2) = 0._WP
         H(3,3) = -1._WP
         H(3,4) = 0._WP
         H(3,5) = 0._WP
         H(3,6) = 0._WP

         H(4,1) = 0._WP
         H(4,2) = 0._WP
         H(4,3) = 1._WP - U
         H(4,4) = -1._WP
         H(4,5) = 0._WP
         H(4,6) = 0._WP

         H(5,1) = 0._WP
         H(5,2) = 0._WP
         H(5,3) = 0._WP
         H(5,4) = 0._WP
         H(5,5) = 1._WP
         H(5,6) = 0._WP

         H(6,1) = 0._WP
         H(6,2) = 0._WP
         H(6,3) = 0._WP
         H(6,4) = 0._WP
         H(6,5) = 0._WP
         H(6,6) = 1._WP

      else

         H(1,1) = 1._WP
         H(1,2) = 0._WP
         H(1,3) = 0._WP
         H(1,4) = 0._WP
         H(1,5) = 0._WP
         H(1,6) = 0._WP

         H(2,1) = 0._WP
         H(2,2) = c_1*omega_c**2
         H(2,3) = 1._WP
         H(2,4) = 0._WP
         H(2,5) = 0._WP
         H(2,6) = 0._WP

         H(3,1) = 0._WP
         H(3,2) = 0._WP
         H(3,3) = -1._WP
         H(3,4) = 0._WP
         H(3,5) = 0._WP
         H(3,6) = 0._WP

         H(4,1) = 0._WP
         H(4,2) = 0._WP
         H(4,3) = 1._WP - U
         H(4,4) = -1._WP
         H(4,5) = 0._WP
         H(4,6) = 0._WP

         H(5,1) = 0._WP
         H(5,2) = 0._WP
         H(5,3) = 0._WP
         H(5,4) = 0._WP
         H(5,5) = 1._WP
         H(5,6) = 0._WP

         H(6,1) = 0._WP
         H(6,2) = 0._WP
         H(6,3) = 0._WP
         H(6,4) = 0._WP
         H(6,5) = 0._WP
         H(6,6) = 1._WP

      end if

    end associate

    ! Finish

    return

  end function H_jcd_

  !****

  function H_lagp_ (this, i, st) result (H)

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: H(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert LAGP variables
    ! to the canonical form

    associate ( &
         V_2 => this%coeff(i,J_V_2))

      ! Set up the matrix
      
      H(1,1) = 1._WP
      H(1,2) = 0._WP
      H(1,3) = 0._WP
      H(1,4) = 0._WP
      H(1,5) = 0._WP
      H(1,6) = 0._WP

      H(2,1) = 1._WP
      H(2,2) = 1._WP/V_2
      H(2,3) = 0._WP
      H(2,4) = 0._WP
      H(2,5) = 0._WP
      H(2,6) = 0._WP

      H(3,1) = 0._WP
      H(3,2) = 0._WP
      H(3,3) = 1._WP
      H(3,4) = 0._WP
      H(3,5) = 0._WP
      H(3,6) = 0._WP

      H(4,1) = 0._WP
      H(4,2) = 0._WP
      H(4,3) = 0._WP
      H(4,4) = 1._WP
      H(4,5) = 0._WP
      H(4,6) = 0._WP

      H(5,1) = 0._WP
      H(5,2) = 0._WP
      H(5,3) = 0._WP
      H(5,4) = 0._WP
      H(5,5) = 1._WP
      H(5,6) = 0._WP

      H(6,1) = 0._WP
      H(6,2) = 0._WP
      H(6,3) = 0._WP
      H(6,4) = 0._WP
      H(6,5) = 0._WP
      H(6,6) = 1._WP

    end associate

    ! Finish

    return

  end function H_lagp_

  !****

  function dH_ (this, i, st) result (dH)

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: dH(this%n_e,this%n_e)

    ! Evaluate the derivative x dH/dx of the transformation matrix H

    select case (this%set)
    case (DZIEM_SET)
       dH = 0._WP
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

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: dH(this%n_e,this%n_e)

    complex(WP) :: omega_c
    complex(WP) :: lambda

    ! Evaluate the derivative x dH/dx of the JCD-variables
    ! transformation matrix H

    associate ( &
         c_1 => this%coeff(i,J_C_1), &
         dc_1 => this%coeff(i,J_DC_1), &
         U => this%coeff(i,J_U), &
         dU => this%coeff(i,J_DU), &
         Omega_rot => this%coeff(i,J_OMEGA_ROT))

      omega_c = this%cx%omega_c(Omega_rot, st)

      lambda = this%cx%lambda(Omega_rot, st)

      ! Set up the matrix (nb: the derivatives of omega_c and lambda is
      ! neglected; this is incorrect when rotation is non-zero)

      ! Set up the matrix

      if (this%l /= 0) then

         dH(1,1) = 0._WP
         dH(1,2) = 0._WP
         dH(1,3) = 0._WP
         dH(1,4) = 0._WP
         dH(1,5) = 0._WP
         dH(1,6) = 0._WP

         dH(2,1) = 0._WP
         dH(2,2) = c_1*dc_1*omega_c**2/lambda
         dH(2,3) = 0._WP
         dH(2,4) = 0._WP
         dH(2,5) = 0._WP
         dH(2,6) = 0._WP
       
         dH(3,1) = 0._WP
         dH(3,2) = 0._WP
         dH(3,3) = 0._WP
         dH(3,4) = 0._WP
         dH(3,5) = 0._WP
         dH(3,6) = 0._WP

         dH(4,1) = 0._WP
         dH(4,2) = 0._WP
         dH(4,3) = -U*dU
         dH(4,4) = 0._WP
         dH(4,5) = 0._WP
         dH(4,6) = 0._WP
         
         dH(5,1) = 0._WP
         dH(5,2) = 0._WP
         dH(5,3) = 0._WP
         dH(5,4) = 0._WP
         dH(5,5) = 0._WP
         dH(5,6) = 0._WP

         dH(6,1) = 0._WP
         dH(6,2) = 0._WP
         dH(6,3) = 0._WP
         dH(6,4) = 0._WP
         dH(6,5) = 0._WP
         dH(6,6) = 0._WP

      else

         dH(1,1) = 0._WP
         dH(1,2) = 0._WP
         dH(1,3) = 0._WP
         dH(1,4) = 0._WP
         dH(1,5) = 0._WP
         dH(1,6) = 0._WP

         dH(2,1) = 0._WP
         dH(2,2) = c_1*dc_1*omega_c**2
         dH(2,3) = 0._WP
         dH(2,4) = 0._WP
         dH(2,5) = 0._WP
         dH(2,6) = 0._WP

         dH(3,1) = 0._WP
         dH(3,2) = 0._WP
         dH(3,3) = 0._WP
         dH(3,4) = 0._WP
         dH(3,5) = 0._WP
         dH(3,6) = 0._WP

         dH(4,1) = 0._WP
         dH(4,2) = 0._WP
         dH(4,3) = U*dU
         dH(4,4) = 0._WP
         dH(4,5) = 0._WP
         dH(4,6) = 0._WP
         
         dH(5,1) = 0._WP
         dH(5,2) = 0._WP
         dH(5,3) = 0._WP
         dH(5,4) = 0._WP
         dH(5,5) = 0._WP
         dH(5,6) = 0._WP

         dH(6,1) = 0._WP
         dH(6,2) = 0._WP
         dH(6,3) = 0._WP
         dH(6,4) = 0._WP
         dH(6,5) = 0._WP
         dH(6,6) = 0._WP

      end if

    end associate

    ! Finish

    return

  end function dH_jcd_

  !****

  function dH_lagp_ (this, i, st) result (dH)

    class(nad_trans_t), intent(in) :: this
    integer, intent(in)            :: i
    class(c_state_t), intent(in)   :: st
    complex(WP)                    :: dH(this%n_e,this%n_e)

    ! Evaluate the derivative x dH/dx of the LAGP-variables
    ! transformation matrix T

    ! Calculate coefficients

    associate ( &
         V_2 => this%coeff(i,J_V_2), &
         dV_2 => this%coeff(i,J_DV_2))

      ! Set up the matrix

      dH(1,1) = 0._WP
      dH(1,2) = 0._WP
      dH(1,3) = 0._WP
      dH(1,4) = 0._WP
      dH(1,5) = 0._WP
      dH(1,6) = 0._WP

      dH(2,1) = 0._WP
      dH(2,2) = -dV_2/V_2
      dH(2,3) = 0._WP
      dH(2,4) = 0._WP
      dH(2,5) = 0._WP
      dH(2,6) = 0._WP

      dH(3,1) = 0._WP
      dH(3,2) = 0._WP
      dH(3,3) = 0._WP
      dH(3,4) = 0._WP
      dH(3,5) = 0._WP
      dH(3,6) = 0._WP

      dH(4,1) = 0._WP
      dH(4,2) = 0._WP
      dH(4,3) = 0._WP
      dH(4,4) = 0._WP
      dH(4,5) = 0._WP
      dH(4,6) = 0._WP

      dH(5,1) = 0._WP
      dH(5,2) = 0._WP
      dH(5,3) = 0._WP
      dH(5,4) = 0._WP
      dH(5,5) = 0._WP
      dH(5,6) = 0._WP

      dH(6,1) = 0._WP
      dH(6,2) = 0._WP
      dH(6,3) = 0._WP
      dH(6,4) = 0._WP
      dH(6,5) = 0._WP
      dH(6,6) = 0._WP

    end associate

    ! Finish

    return

  end function dH_lagp_

end module gyre_nad_trans
