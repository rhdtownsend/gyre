! Module   : gyre_ad_bvp
! Purpose  : solve adiabatic BVPs
!
! Copyright 2013 Rich Townsend
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
  use gyre_mech_coeffs
  use gyre_oscpar
  use gyre_gridpar
  use gyre_numpar
  use gyre_ad_shooter
  use gyre_ad_bound
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(bvp_t) :: ad_bvp_t
     private
     class(mech_coeffs_t), allocatable, public :: mc
     type(oscpar_t), public                    :: op
     type(gridpar_t), public                   :: gp
     type(numpar_t), public                    :: np
     type(ad_shooter_t)                        :: sh
     type(ad_bound_t)                          :: bd
     type(sysmtx_t)                            :: sm
     real(WP), allocatable                     :: x(:)
     integer                                   :: e_norm
     integer, public                           :: n
     integer, public                           :: n_e
   contains 
     private
     procedure, public :: init
     procedure, public :: set_norm
     procedure, public :: discrim
     procedure         :: build
     procedure, public :: recon
  end type ad_bvp_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_bp
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: ad_bvp_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init (this, mc, op, gp, np, x)

    class(ad_bvp_t), intent(out)     :: this
    class(mech_coeffs_t), intent(in) :: mc
    type(oscpar_t), intent(in)       :: op
    type(gridpar_t), intent(in)      :: gp
    type(numpar_t), intent(in)       :: np
    real(WP), intent(in)             :: x(:)

    integer :: n

    ! Initialize the ad_bvp

    allocate(this%mc, SOURCE=mc)
    this%op = op
    this%gp = gp
    this%np = np

    call this%sh%init(this%mc, this%op, this%np)
    call this%bd%init(this%mc, this%op)

    n = SIZE(x)

    call this%sm%init(n-1, this%sh%n_e, this%bd%n_i, this%bd%n_o)

    this%x = x

    this%e_norm = 0

    this%n = n
    this%n_e = this%sh%n_e

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_bp (this, root_rank)

    class(ad_bvp_t), intent(inout) :: this
    integer, intent(in)            :: root_rank

    class(mech_coeffs_t), allocatable :: mc
    type(oscpar_t)                    :: op
    type(gridpar_t)                   :: gp
    type(numpar_t)                    :: np
    real(WP), allocatable             :: x(:)

    ! Broadcast the bvp

    if(MPI_RANK == root_rank) then

       call bcast_alloc(this%mc, root_rank)
      
       call bcast(this%op, root_rank)
       call bcast(this%gp, root_rank)
       call bcast(this%np, root_rank)

       call bcast_alloc(this%x, root_rank)

    else

       call bcast_alloc(mc, root_rank)

       call bcast(op, root_rank)
       call bcast(gp, root_rank)
       call bcast(np, root_rank)
    
       call bcast_alloc(x, root_rank)

       call this%init(mc, op, gp, np, x)

    endif

    ! Finish

    return

  end subroutine bcast_bp

  $endif

!****

  subroutine set_norm (this, omega)

    class(ad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)        :: omega

    type(ext_complex_t) :: discrim

    ! Evaluate the discriminant

    discrim = this%discrim(omega)

    ! Set the normalizing exponent based on this discriminant

    this%e_norm = discrim%e

    ! Finish

    return

  end subroutine set_norm

!****

  function discrim (this, omega, norm)

    class(ad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)        :: omega
    logical, intent(in), optional  :: norm
    type(ext_complex_t)            :: discrim

    logical :: norm_

    if(PRESENT(norm)) then
       norm_ = norm
    else
       norm_ = .FALSE.
    endif

    ! Evaluate the discriminant as the determinant of the sysmtx

    call this%build(omega)

    discrim = this%sm%determinant()

    ! Apply the normalization

    if(norm_) discrim%e = discrim%e - this%e_norm

    ! Finish

    return

  end function discrim

!****

  subroutine build (this, omega)

    class(ad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)        :: omega

    ! Set up the sysmtx

    call this%sm%set_inner_bound(this%bd%inner_bound(omega))
    call this%sm%set_outer_bound(this%bd%outer_bound(omega))

    call this%sh%shoot(omega, this%x, this%sm)

    ! Finish

    return

  end subroutine build

!****

  subroutine recon (this, omega, x, y)

    class(ad_bvp_t), intent(inout)        :: this
    complex(WP), intent(in)               :: omega
    real(WP), allocatable, intent(out)    :: x(:)
    complex(WP), allocatable, intent(out) :: y(:,:)

    complex(WP) :: y_sh(this%n_e,this%n)
    integer     :: dn(this%n-1)

    ! Reconstruct the solution on the shooting grid

    call this%build(omega)

    y_sh = RESHAPE(this%sm%null_vector(), SHAPE(y_sh))

    ! Set up the grid

    select case (this%gp%grid_type)
    case ('GEOM')
       call build_geom_grid(this%gp%s, this%gp%n_grid, x)
    case ('LOG')
       call build_geom_grid(this%gp%s, this%gp%n_grid, x)
    case ('CLONE')
       x = this%x
    case ('DISP')
       dn = 0
       call plan_dispersion_grid(this%x, this%mc, omega, this%op, &
                                 this%gp%alpha_osc, this%gp%alpha_exp, this%gp%n_center, this%gp%n_floor, dn)
       call build_oversamp_grid(this%x, dn, x)
    case default
       $ABORT(Invalid grid_type)
    end select

    ! Reconstruct the full solution

    allocate(y(this%n_e,SIZE(x)))

    call this%sh%recon(omega, this%x, y_sh, x, y)

    ! Finish

    return

  end subroutine recon

end module gyre_ad_bvp
