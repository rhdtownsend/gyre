! Module   : gyre_ad_bvp
! Purpose  : boundary-value solver (adiabatic)
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

module gyre_ad_bvp

  ! Uses

  use core_kinds

  use gyre_ad_bound
  use gyre_ad_eqns
  use gyre_ad_vars
  use gyre_bvp
  use gyre_ext
  use gyre_ivp
  use gyre_ivp_factory
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_num_par
  use gyre_osc_par
  use gyre_sysmtx
  use gyre_sysmtx_factory
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_bvp_t) :: ad_bvp_t
     class(model_t), pointer :: ml => null()
     type(ad_vars_t)         :: vr
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

  function ad_bvp_t_ (x, ml, mp, op, np, omega_min, omega_max) result (bp)

    real(WP), intent(in)                :: x(:)
    class(model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)        :: mp
    type(osc_par_t), intent(in)         :: op
    type(num_par_t), intent(in)         :: np
    real(WP), intent(in)                :: omega_min
    real(WP), intent(in)                :: omega_max
    type(ad_bvp_t)                      :: bp

    class(r_rot_t), allocatable    :: rt
    type(ad_eqns_t)                :: eq
    integer                        :: n
    real(WP)                       :: x_i
    real(WP)                       :: x_o
    type(ad_bound_t)               :: bd
    class(r_ivp_t), allocatable    :: iv
    class(r_sysmtx_t), allocatable :: sm

    ! Construct the ad_bvp_t

    ! Initialize the rotational effects

    allocate(rt, SOURCE=r_rot_t(ml, mp, op))
 
    ! Initialize the equations

    eq = ad_eqns_t(ml, rt, op)

    ! Initialize the boundary conditions

    n = SIZE(x)

    x_i = x(1)
    x_o = x(n)

    bd = ad_bound_t(ml, rt, op, x_i, x_o)

    ! Initialize the IVP solver

    allocate(iv, SOURCE=r_ivp_t(eq, np))

    ! Initialize the system matrix

    allocate(sm, SOURCE=r_sysmtx_t(n-1, eq%n_e, bd%n_i, bd%n_o, np))

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t(x, bd, iv, sm, omega_min, omega_max)

    bp%ml => ml
    bp%vr = ad_vars_t(ml, rt, op)
    
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
       y(1:4,i) = MATMUL(this%vr%T(x(i), omega), y_(:,i))
       y(5:6,i) = 0._WP
    end do

    y_ref(1:4) = MATMUL(this%vr%T(x_ref, omega), y_ref_)
    y_ref(5:6) = 0._WP

    ! Finish

    return

  end subroutine recon_

end module gyre_ad_bvp

