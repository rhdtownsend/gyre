! Module   : gyre_ad_bvp
! Purpose  : boundary-value solver (adiabatic)
!
! Copyright 2013-2014 Rich Townsend
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

module gyre_ad_bvp

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_ext
  use gyre_ivp
  use gyre_mode
  use gyre_modepar
  use gyre_model
  use gyre_numpar
  use gyre_oscpar
  use gyre_sysmtx
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_bvp_t) :: ad_bvp_t
   contains
     private
     procedure, public :: recon => recon_
  end type ad_bvp_t

  ! Interfaces

  interface ad_bvp_t
     module procedure ad_bvp_t_
  end interface ad_bvp_t

  ! Access specifiers

  private

  public :: ad_bvp_t

  ! Procedures

contains

  function ad_bvp_t_ (x, ml, mp, op, np) result (bp)

    use gyre_dopp_rot
    use gyre_null_rot
    use gyre_trad_rot

    use gyre_ad_jacob
    use gyre_ad_bound

    use gyre_magnus_ivp
    use gyre_colloc_ivp
    use gyre_findiff_ivp

    use gyre_block_sysmtx

    real(WP), intent(in)                :: x(:)
    class(model_t), pointer, intent(in) :: ml
    type(modepar_t), intent(in)         :: mp
    type(oscpar_t), intent(in)          :: op
    type(numpar_t), intent(in)          :: np
    type(ad_bvp_t), target              :: bp

    class(r_rot_t), allocatable    :: rt
    type(ad_jacob_t)               :: jc
    integer                        :: n
    real(WP)                       :: x_i
    real(WP)                       :: x_o
    type(ad_bound_t)               :: bd
    class(r_ivp_t), allocatable    :: iv
    class(r_sysmtx_t), allocatable :: sm

    ! Construct the ad_bvp_t

    ! Initialize the rotational effects

    select case (op%rotation_method)
    case ('DOPPLER')
       allocate(rt, SOURCE=r_dopp_rot_t(ml, mp))
    case ('NULL')
       allocate(rt, SOURCE=r_null_rot_t(mp))
    case ('TRAD')
       allocate(rt, SOURCE=r_trad_rot_t(ml, mp))
    case default
       $ABORT(Invalid rotation_method)
    end select
 
    ! Initialize the jacobian

    jc = ad_jacob_t(ml, rt, op%variables_type)

    ! Initialize the boundary conditions

    n = SIZE(x)

    x_i = x(1)
    x_o = x(n)

    bd = ad_bound_t(ml, rt, jc, x_i, x_o, op%inner_bound_type, op%outer_bound_type)

    ! Initialize the IVP solver

    select case (np%ivp_solver_type)
    case ('MAGNUS_GL2')
       allocate(iv, SOURCE=r_magnus_ivp_t(jc, 'GL2'))
    case ('MAGNUS_GL4')
       allocate(iv, SOURCE=r_magnus_ivp_t(jc, 'GL4'))
    case ('MAGNUS_GL6')
       allocate(iv, SOURCE=r_magnus_ivp_t(jc, 'GL6'))
    case ('COLLOC_GL2')
       allocate(iv, SOURCE=r_colloc_ivp_t(jc, 'GL2'))
    case ('COLLOC_GL4')
       allocate(iv, SOURCE=r_colloc_ivp_t(jc, 'GL4'))
    case ('FINDIFF')
       allocate(iv, SOURCE=r_findiff_ivp_t(jc))
    case default
       $ABORT(Invalid ivp_solver_type)
    end select

    ! Initialize the system matrix

    if (np%use_banded) then
       $ABORT(Not yet implemented)
    else
       allocate(sm, SOURCE=r_block_sysmtx_t(n-1, jc%n_e, bd%n_i, bd%n_o))
    endif

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t(x, ml, jc, bd, iv, sm)

    ! Finish

    return

  end function ad_bvp_t_

!****

  subroutine recon_ (this, omega, x, x_ref, y, y_ref, discrim)

    class(ad_bvp_t), intent(inout) :: this
    real(WP), intent(in)           :: omega
    real(WP), intent(in)           :: x(:)
    real(WP), intent(in)           :: x_ref
    real(WP), intent(out)          :: y(:,:)
    real(WP), intent(out)          :: y_ref(:)
    type(r_ext_t), intent(out)     :: discrim

    real(WP) :: y_(4,SIZE(x))
    real(WP) :: y_ref_(4)
    integer  :: n
    integer  :: i

    $CHECK_BOUNDS(SIZE(y, 1),6)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    $CHECK_BOUNDS(SIZE(y_ref),6)

    ! Reconstruct the solution

    call this%r_bvp_t%recon(omega, x, x_ref, y_, y_ref_, discrim)

    ! Convert to the canonical (6-variable) solution

    n = SIZE(x)

    !$OMP PARALLEL DO 
    do i = 1, n
       y(1:4,i) = MATMUL(this%jc%T(x(i), omega, .TRUE.), y_(:,i))
       y(5:6,i) = 0._WP
    end do

    y_ref(1:4) = MATMUL(this%jc%T(x_ref, omega, .TRUE.), y_ref_)
    y_ref(5:6) = 0._WP

    ! Finish

    return

  end subroutine recon_

end module gyre_ad_bvp

