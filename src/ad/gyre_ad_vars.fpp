! Module   : gyre_ad_vars
! Purpose  : adiabatic variables transformations
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

module gyre_ad_vars

  ! Uses

  use core_kinds

  use gyre_grid
  use gyre_linalg
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: GYRE_SET = 0
  integer, parameter :: DZIEM_SET = 1
  integer, parameter :: JCD_SET = 2
  integer, parameter :: MIX_SET = 3
  integer, parameter :: LAGP_SET = 4

  ! Derived-type definitions

  type :: ad_vars_t
     private
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     integer                     :: l
     integer                     :: set
   contains
     private
     procedure, public :: G
     procedure         :: G_dziem_
     procedure         :: G_jcd_
     procedure         :: G_mix_
     procedure         :: G_lagp_
     procedure, public :: H
     procedure         :: H_dziem_
     procedure         :: H_jcd_
     procedure         :: H_mix_
     procedure         :: H_lagp_
     procedure, public :: dH
     procedure         :: dH_jcd_
     procedure         :: dH_mix_
     procedure         :: dH_lagp_
  end type ad_vars_t

  ! Interfaces

  interface ad_vars_t
     module procedure ad_vars_t_
  end interface ad_vars_t

  ! Access specifiers

  private

  public :: ad_vars_t

  ! Procedures

contains

  function ad_vars_t_ (ml, gr, md_p, os_p) result (vr)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(ad_vars_t)                     :: vr

    ! Construct the ad_vars_t

    call check_model(ml, [I_V_2,I_AS,I_U,I_C_1,I_GAMMA_1])

    vr%ml => ml

    allocate(vr%rt, SOURCE=r_rot_t(ml, gr, md_p, os_p))

    select case (os_p%variables_set)
    case ('GYRE')
       vr%set = GYRE_SET
    case ('DZIEM')
       vr%set = DZIEM_SET
    case ('JCD')
       vr%set = JCD_SET
    case ('MIX')
       vr%set = MIX_SET
    case ('LAGP')
       vr%set = LAGP_SET
    case default
       $ABORT(Invalid variables_set)
    end select

    vr%l = md_p%l

    ! Finish

    return

  end function ad_vars_t_

  !****

  function G (this, pt, omega)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: G(4,4)

    ! Evaluate the transformation matrix to convert variables from
    ! the canonical form

    select case (this%set)
    case (GYRE_SET)
       G = identity_matrix(4)
    case (DZIEM_SET)
       G = this%G_dziem_(pt, omega)
    case (JCD_SET)
       G = this%G_jcd_(pt, omega)
    case (MIX_SET)
       G = this%G_mix_(pt, omega)
    case (LAGP_SET)
       G = this%G_lagp_(pt, omega)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function G

  !****

  function G_dziem_ (this, pt, omega) result (G)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: G(4,4)

    ! Evaluate the transformation matrix to convert DZIEM variables
    ! from the canonical form

    ! Set up the matrix
      
    G(1,1) = 1._WP
    G(1,2) = 0._WP
    G(1,3) = 0._WP
    G(1,4) = 0._WP
          
    G(2,1) = 0._WP
    G(2,2) = 1._WP
    G(2,3) = 1._WP
    G(2,4) = 0._WP
          
    G(3,1) = 0._WP
    G(3,2) = 0._WP
    G(3,3) = 1._WP
    G(3,4) = 0._WP

    G(4,1) = 0._WP
    G(4,2) = 0._WP
    G(4,3) = 0._WP
    G(4,4) = 1._WP

    ! Finish

    return

  end function G_dziem_

  !****

  function G_jcd_ (this, pt, omega) result (G)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: G(4,4)

    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! from the canonical form

    ! Calculate coefficients

    U = this%ml%coeff(I_U, pt)
    c_1 = this%ml%coeff(I_C_1, pt)

    lambda = this%rt%lambda(pt, omega)

    omega_c = this%rt%omega_c(pt, omega)

    ! Set up the matrix
      
    if (this%l /= 0) then

       G(1,1) = 1._WP
       G(1,2) = 0._WP
       G(1,3) = 0._WP
       G(1,4) = 0._WP
          
       G(2,1) = 0._WP
       G(2,2) = lambda/(c_1*omega_c**2)
       G(2,3) = lambda/(c_1*omega_c**2)
       G(2,4) = 0._WP
          
       G(3,1) = 0._WP
       G(3,2) = 0._WP
       G(3,3) = -1._WP
       G(3,4) = 0._WP

       G(4,1) = 0._WP
       G(4,2) = 0._WP
       G(4,3) = -(1._WP - U)
       G(4,4) = -1._WP

    else

       G(1,1) = 1._WP
       G(1,2) = 0._WP
       G(1,3) = 0._WP
       G(1,4) = 0._WP

       G(2,1) = 0._WP
       G(2,2) = 1._WP/(c_1*omega_c**2)
       G(2,3) = 1._WP/(c_1*omega_c**2)
       G(2,4) = 0._WP

       G(3,1) = 0._WP
       G(3,2) = 0._WP
       G(3,3) = -1._WP
       G(3,4) = 0._WP

       G(4,1) = 0._WP
       G(4,2) = 0._WP
       G(4,3) = -(1._WP - U)
       G(4,4) = -1._WP

    endif

    ! Finish

    return

  end function G_jcd_

  !****

  function G_mix_ (this, pt, omega) result (G)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: G(4,4)

    real(WP) :: U

    ! Evaluate the transformation matrix to convert MIX variables
    ! from the canonical form

    ! Calculate coefficients

    U = this%ml%coeff(I_U, pt)

    ! Set up the matrix

    G(1,1) = 1._WP
    G(1,2) = 0._WP
    G(1,3) = 0._WP
    G(1,4) = 0._WP

    G(2,1) = 0._WP
    G(2,2) = 1._WP
    G(2,3) = 0._WP
    G(2,4) = 0._WP

    G(3,1) = 0._WP
    G(3,2) = 0._WP
    G(3,3) = -1._WP
    G(3,4) = 0._WP

    G(4,1) = 0._WP
    G(4,2) = 0._WP
    G(4,3) = -(1._WP - U)
    G(4,4) = -1._WP

    ! Finish

    return

  end function G_mix_

  !****

  function G_lagp_ (this, pt, omega) result (G)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: G(4,4)

    real(WP) :: V_2

    $ASSERT(.NOT. this%ml%is_vacuum(pt),Cannot use LAGP variables at vacuum points)

    ! Evaluate the transformation matrix to convert LAGP variables
    ! from the canonical form

    ! Calculate coefficients

    V_2 = this%ml%coeff(I_V_2, pt)

    ! Set up the matrix

    G(1,1) = 1._WP
    G(1,2) = 0._WP
    G(1,3) = 0._WP
    G(1,4) = 0._WP

    G(2,1) = -V_2
    G(2,2) = V_2
    G(2,3) = 0._WP
    G(2,4) = 0._WP

    G(3,1) = 0._WP
    G(3,2) = 0._WP
    G(3,3) = 1._WP
    G(3,4) = 0._WP

    G(4,1) = 0._WP
    G(4,2) = 0._WP
    G(4,3) = 0._WP
    G(4,4) = 1._WP

    ! Finish

    return

  end function G_lagp_

  !****

  function H (this, pt, omega)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: H(4,4)

    ! Evaluate the transformation matrix to convert variables to
    ! canonical form

    select case (this%set)
    case (GYRE_SET)
       H = identity_matrix(4)
    case (DZIEM_SET)
       H = this%H_dziem_(pt, omega)
    case (JCD_SET)
       H = this%H_jcd_(pt, omega)
    case (MIX_SET)
       H = this%H_mix_(pt, omega)
    case (LAGP_SET)
       H = this%H_lagp_(pt, omega)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function H

  !****

  function H_dziem_ (this, pt, omega) result (H)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: H(4,4)

    ! Evaluate the transformation matrix to convert DZIEM variables
    ! to the canonical form

    ! Set up the matrix
      
    H(1,1) = 1._WP
    H(1,2) = 0._WP
    H(1,3) = 0._WP
    H(1,4) = 0._WP
       
    H(2,1) = 0._WP
    H(2,2) = 1._WP
    H(2,3) = -1._WP
    H(2,4) = 0._WP

    H(3,1) = 0._WP
    H(3,2) = 0._WP
    H(3,3) = 1._WP
    H(3,4) = 0._WP

    H(4,1) = 0._WP
    H(4,2) = 0._WP
    H(4,3) = 0._WP
    H(4,4) = 1._WP
    
    ! Finish

    return

  end function H_dziem_

  !****

  function H_jcd_ (this, pt, omega) result (H)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: H(4,4)

    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! to the canonical form

    ! Calculate coefficients

    U = this%ml%coeff(I_U, pt)
    c_1 = this%ml%coeff(I_C_1, pt)

    lambda = this%rt%lambda(pt, omega)

    omega_c = this%rt%omega_c(pt, omega)

    ! Set up the matrix
      
    if (this%l /= 0._WP) then

       H(1,1) = 1._WP
       H(1,2) = 0._WP
       H(1,3) = 0._WP
       H(1,4) = 0._WP
       
       H(2,1) = 0._WP
       H(2,2) = c_1*omega_c**2/lambda
       H(2,3) = 1._WP
       H(2,4) = 0._WP

       H(3,1) = 0._WP
       H(3,2) = 0._WP
       H(3,3) = -1._WP
       H(3,4) = 0._WP

       H(4,1) = 0._WP
       H(4,2) = 0._WP
       H(4,3) = 1._WP - U
       H(4,4) = -1._WP

    else

       H(1,1) = 1._WP
       H(1,2) = 0._WP
       H(1,3) = 0._WP
       H(1,4) = 0._WP

       H(2,1) = 0._WP
       H(2,2) = c_1*omega_c**2
       H(2,3) = 1._WP
       H(2,4) = 0._WP

       H(3,1) = 0._WP
       H(3,2) = 0._WP
       H(3,3) = -1._WP
       H(3,4) = 0._WP

       H(4,1) = 0._WP
       H(4,2) = 0._WP
       H(4,3) = 1._WP - U
       H(4,4) = -1._WP

    endif

    ! Finish

    return

  end function H_jcd_

  !****

  function H_mix_ (this, pt, omega) result (H)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: H(4,4)

    real(WP) :: U

    ! Evaluate the transformation matrix to convert MIX variables
    ! to the canonical form

    ! Calculate coefficients

    U = this%ml%coeff(I_U, pt)

    ! Set up the matrix

    H(1,1) = 1._WP
    H(1,2) = 0._WP
    H(1,3) = 0._WP
    H(1,4) = 0._WP

    H(2,1) = 0._WP
    H(2,2) = 1._WP
    H(2,3) = 0._WP
    H(2,4) = 0._WP

    H(3,1) = 0._WP
    H(3,2) = 0._WP
    H(3,3) = -1._WP
    H(3,4) = 0._WP

    H(4,1) = 0._WP
    H(4,2) = 0._WP
    H(4,3) = 1._WP - U
    H(4,4) = -1._WP

    ! Finish

    return

  end function H_mix_

  !****

  function H_lagp_ (this, pt, omega) result (H)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: H(4,4)

    real(WP) :: V_2

    $ASSERT(.NOT. this%ml%is_vacuum(pt),Cannot use LAGP variables at vacuum points)

    ! Evaluate the transformation matrix to convert LAGP variables
    ! to the canonical form

    ! Calculate coefficients

    V_2 = this%ml%coeff(I_V_2, pt)

    ! Set up the matrix

    H(1,1) = 1._WP
    H(1,2) = 0._WP
    H(1,3) = 0._WP
    H(1,4) = 0._WP

    H(2,1) = 1._WP
    H(2,2) = 1._WP/V_2
    H(2,3) = 0._WP
    H(2,4) = 0._WP

    H(3,1) = 0._WP
    H(3,2) = 0._WP
    H(3,3) = 1._WP
    H(3,4) = 0._WP

    H(4,1) = 0._WP
    H(4,2) = 0._WP
    H(4,3) = 0._WP
    H(4,4) = 1._WP

    ! Finish

    return

  end function H_lagp_

  !****

  function dH (this, pt, omega)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: dH(4,4)

    ! Evaluate the derivative x dH/dx of the transformation matrix H

    select case (this%set)
    case (GYRE_SET)
       dH = 0._WP
    case (DZIEM_SET)
       dH = 0._WP
    case (JCD_SET)
       dH = this%dH_jcd_(pt, omega)
    case (MIX_SET)
       dH = this%dH_mix_(pt, omega)
    case (LAGP_SET)
       dH = this%dH_lagp_(pt, omega)
    case default
       $ABORT(Invalid set)
    end select

    ! Finish

    return

  end function dH

  !****

  function dH_jcd_ (this, pt, omega) result (dH)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: dH(4,4)

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: omega_c

    ! Evaluate the derivative x dH/dx of the JCD-variables
    ! transformation matrix H

    ! Calculate coefficients

    V_g = this%ml%coeff(I_V_2, pt)*pt%x**2/this%ml%coeff(I_GAMMA_1, pt)
    As = this%ml%coeff(I_AS, pt)
    U = this%ml%coeff(I_U, pt)
    c_1 = this%ml%coeff(I_C_1, pt)

    lambda = this%rt%lambda(pt, omega)

    omega_c = this%rt%omega_c(pt, omega)

    ! Set up the matrix (nb: the derivatives of omega_c and lambda are
    ! neglected; this is incorrect when rotation is non-zero)
      
    if (this%l /= 0._WP) then

       dH(1,1) = 0._WP
       dH(1,2) = 0._WP
       dH(1,3) = 0._WP
       dH(1,4) = 0._WP
       
       dH(2,1) = 0._WP
       dH(2,2) = c_1*(3._WP - U)*omega_c**2/lambda
       dH(2,3) = 0._WP
       dH(2,4) = 0._WP

       dH(3,1) = 0._WP
       dH(3,2) = 0._WP
       dH(3,3) = 0._WP
       dH(3,4) = 0._WP

       dH(4,1) = 0._WP
       dH(4,2) = 0._WP
       dH(4,3) = U*(V_g + As + U - 3._WP)
       dH(4,4) = 0._WP

    else

       dH(1,1) = 0._WP
       dH(1,2) = 0._WP
       dH(1,3) = 0._WP
       dH(1,4) = 0._WP

       dH(2,1) = 0._WP
       dH(2,2) = c_1*(3._WP - U)*omega_c**2
       dH(2,3) = 0._WP
       dH(2,4) = 0._WP

       dH(3,1) = 0._WP
       dH(3,2) = 0._WP
       dH(3,3) = 0._WP
       dH(3,4) = 0._WP

       dH(4,1) = 0._WP
       dH(4,2) = 0._WP
       dH(4,3) = U*(V_g + As + U - 3._WP)
       dH(4,4) = 0._WP

    endif

    ! Finish

    return

  end function dH_jcd_

  !****

  function dH_mix_ (this, pt, omega) result (dH)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: dH(4,4)

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U

    ! Evaluate the derivative x dH/dx of the MIX-variables
    ! transformation matrix H

    ! Calculate coefficients

    V_g = this%ml%coeff(I_V_2, pt)*pt%x**2/this%ml%coeff(I_GAMMA_1, pt)
    As = this%ml%coeff(I_AS, pt)
    U = this%ml%coeff(I_U, pt)

    ! Set up the matrix

    dH(1,1) = 0._WP
    dH(1,2) = 0._WP
    dH(1,3) = 0._WP
    dH(1,4) = 0._WP

    dH(2,1) = 0._WP
    dH(2,2) = 0._WP
    dH(2,3) = 0._WP
    dH(2,4) = 0._WP

    dH(3,1) = 0._WP
    dH(3,2) = 0._WP
    dH(3,3) = 0._WP
    dH(3,4) = 0._WP

    dH(4,1) = 0._WP
    dH(4,2) = 0._WP
    dH(4,3) = U*(V_g + As + U - 3._WP)
    dH(4,4) = 0._WP

    ! Finish

    return

  end function dH_mix_

!****

  function dH_lagp_ (this, pt, omega) result (dH)

    class(ad_vars_t), intent(in) :: this
    type(point_t), intent(in)    :: pt
    real(WP), intent(in)         :: omega
    real(WP)                     :: dH(4,4)

    real(WP) :: V_2
    real(WP) :: V
    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U

    $ASSERT(.NOT. this%ml%is_vacuum(pt),Cannot use LAGP variables at vacuum points)

    ! Evaluate the derivative x dH/dx of the LAGP-variables
    ! transformation matrix H

    ! Calculate coefficients

    V_2 = this%ml%coeff(I_V_2, pt)
    V = V_2*pt%x**2
    V_g = V/this%ml%coeff(I_GAMMA_1, pt)
    As = this%ml%coeff(I_AS, pt)
    U = this%ml%coeff(I_U, pt)

    ! Set up the matrix

    dH(1,1) = 0._WP
    dH(1,2) = 0._WP
    dH(1,3) = 0._WP
    dH(1,4) = 0._WP

    dH(2,1) = 0._WP
    dH(2,2) = -(-V_g - As + U + V - 3)/V_2
    dH(2,3) = 0._WP
    dH(2,4) = 0._WP

    dH(3,1) = 0._WP
    dH(3,2) = 0._WP
    dH(3,3) = 0._WP
    dH(3,4) = 0._WP

    dH(4,1) = 0._WP
    dH(4,2) = 0._WP
    dH(4,3) = 0._WP
    dH(4,4) = 0._WP

    ! Finish

    return

  end function dH_lagp_

end module gyre_ad_vars
