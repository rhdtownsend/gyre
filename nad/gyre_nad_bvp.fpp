! Module   : gyre_nad_bvp
! Purpose  : solve nonadiabatic BVPs
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

module gyre_nad_bvp

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_bvp
  use gyre_base_coeffs
  use gyre_therm_coeffs
  $if($MPI)
  use gyre_base_coeffs_mpi
  use gyre_therm_coeffs_mpi
  $endif
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_nad_shooter
  use gyre_nad_bound
  use gyre_nad_jacobian
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_grid
  use gyre_mode

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(bvp_t) :: nad_bvp_t
     private
     class(base_coeffs_t), allocatable  :: bc
     class(therm_coeffs_t), allocatable :: tc
     type(oscpar_t)                     :: op
     type(numpar_t)                     :: np
     type(gridpar_t), allocatable       :: shoot_gp(:)
     type(gridpar_t), allocatable       :: recon_gp(:)
     type(nad_shooter_t)                :: sh
     type(nad_bound_t)                  :: bd
     type(sysmtx_t)                     :: sm
     real(WP), allocatable              :: x_in(:)
     real(WP), allocatable              :: x(:)
     real(WP)                           :: x_ad
     integer                            :: e_norm
     integer, public                    :: n
     integer, public                    :: n_e
   contains 
     private
     procedure, public :: init
!     procedure, public :: get_bc
!     procedure, public :: get_tc
!     procedure, public :: get_op
     procedure, public :: get_np
!     procedure, public :: get_gp
     procedure, public :: set_norm
     procedure, public :: set_x_ad
     procedure, public :: discrim
     procedure         :: build
     procedure         :: recon
     procedure, public :: mode
  end type nad_bvp_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_bp
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: nad_bvp_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init (this, bc, op, np, shoot_gp, recon_gp, x_in, tc)

    class(nad_bvp_t), intent(out)     :: this
    class(base_coeffs_t), intent(in)  :: bc
    type(oscpar_t), intent(in)        :: op
    type(numpar_t), intent(in)        :: np
    type(gridpar_t), intent(in)       :: shoot_gp(:)
    type(gridpar_t), intent(in)       :: recon_gp(:)
    real(WP), allocatable, intent(in) :: x_in(:)
    class(therm_coeffs_t), intent(in) :: tc

    integer               :: n
    real(WP), allocatable :: x_cc(:)

    ! Initialize the nad_bvp

    ! Create the shooting grid

    call build_grid(shoot_gp, bc, op, x_in, this%x, tc)

    n = SIZE(this%x)

    ! Set up components
    
    allocate(this%bc, SOURCE=bc)
    allocate(this%tc, SOURCE=tc)

    this%op = op
    this%np = np

    this%shoot_gp = shoot_gp
    this%recon_gp = recon_gp

    call this%sh%init(this%bc, this%tc, this%op, this%np)
    call this%bd%init(this%bc, this%tc, this%op)

    call this%sm%init(n-1, this%sh%n_e, this%bd%n_i, this%bd%n_o)

    if(ALLOCATED(x_in)) this%x_in = x_in

    this%x_ad = 0._WP

    this%e_norm = 0

    this%n = n
    this%n_e = this%sh%n_e

    ! Set up the coefficient caches

    ! Set up the coefficient caches

    x_cc = [this%x(1),this%sh%abscissa(this%x),this%x(n)]

    call this%bc%fill_cache(x_cc)
    call this%tc%fill_cache(x_cc)

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_bp (this, root_rank)

    class(nad_bvp_t), intent(inout) :: this
    integer, intent(in)             :: root_rank

    class(base_coeffs_t), allocatable  :: bc
    class(therm_coeffs_t), allocatable :: tc
    type(oscpar_t)                     :: op
    type(numpar_t)                     :: np
    type(gridpar_t), allocatable       :: shoot_gp(:)
    type(gridpar_t), allocatable       :: recon_gp(:)
    real(WP), allocatable              :: x_in(:)

    ! Broadcast the bvp

    if(MPI_RANK == root_rank) then

       call bcast_alloc(this%bc, root_rank)
       call bcast_alloc(this%tc, root_rank)
      
       call bcast(this%op, root_rank)
       call bcast(this%np, root_rank)

       call bcast_alloc(this%shoot_gp, root_rank)
       call bcast_alloc(this%recon_gp, root_rank)

       call bcast_alloc(this%x_in, root_rank)

    else

       call bcast_alloc(bc, root_rank)
       call bcast_alloc(tc, root_rank)

       call bcast(op, root_rank)
       call bcast(np, root_rank)

       call bcast_alloc(shoot_gp, root_rank)
       call bcast_alloc(recon_gp, root_rank)

       call bcast_alloc(x_in, root_rank)

       call this%init(bc, tc, op, np, shoot_gp, recon_gp, x_in)

    endif

    ! Finish

    return

  end subroutine bcast_bp

  $endif

!****

  $define $GET $sub

  $local $ITEM_NAME $1
  $local $ITEM_TYPE $2

  function get_$ITEM_NAME (this) result (item)

    class(nad_bvp_t), intent(in), target :: this
    $ITEM_TYPE, pointer                  :: item

    ! Get a pointer to the item

    item => this%$ITEM_NAME

    ! Finish

    return

  end function get_$ITEM_NAME
  
  $endsub

!  $GET(bc,class(base_coeffs_t))
!  $GET(tc,class(therm_coeffs_t))
!  $GET(op,type(oscpar_t))
!  $GET(gp,type(gridpar_t))
  $GET(np,type(numpar_t))
    
!****

  subroutine set_norm (this, omega)

    class(nad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega

    type(ext_complex_t) :: discrim

    ! Evaluate the discriminant

    discrim = this%discrim(omega)

    ! Set the normalizing exponent based on this discriminant

    this%e_norm = exponent(discrim)

    ! Finish

    return

  end subroutine set_norm

!****

  subroutine set_x_ad (this, omega)

    class(nad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega

    integer :: k

    ! Decide where to switch from the adiabatic equations (interior)
    ! to the non-adiabastic ones (exterior)

    this%x_ad = 0._WP

    x_ad_loop : do k = this%n,2,-1

       if(this%tc%tau_thm(this%x(k))*REAL(omega)*this%np%theta_ad > 1._WP) then
          this%x_ad = this%x(k)
          exit x_ad_loop
       endif

    end do x_ad_loop

    ! Finish

    return

  end subroutine set_x_ad

!****

  function discrim (this, omega, norm)

    class(nad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega
    logical, intent(in), optional   :: norm
    type(ext_complex_t)             :: discrim

    logical :: norm_

    if(PRESENT(norm)) then
       norm_ = norm
    else
       norm_ = .FALSE.
    endif

    ! Evaluate the discriminant as the determinant of the sysmtx

    call this%build(omega)

    call this%sm%determinant(discrim)

    ! Apply the normalization

    if(norm_) discrim = scale(discrim, -this%e_norm)

    ! Finish

    return

  end function discrim

!****

  subroutine build (this, omega)

    class(nad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega

    ! Set up the sysmtx

    call this%bc%enable_cache()
    call this%tc%enable_cache()

    call this%sm%set_inner_bound(this%bd%inner_bound(this%x(1), omega))
    call this%sm%set_outer_bound(this%bd%outer_bound(this%x(this%n), omega))

    call this%sh%shoot(omega, this%x, this%sm, this%x_ad)

    call this%bc%disable_cache()
    call this%tc%disable_cache()

    ! Finish

    return

  end subroutine build

!****

  subroutine recon (this, omega, x, y, discrim)

    class(nad_bvp_t), intent(inout)       :: this
    complex(WP), intent(in)               :: omega
    real(WP), allocatable, intent(out)    :: x(:)
    complex(WP), allocatable, intent(out) :: y(:,:)
    type(ext_complex_t), intent(out)      :: discrim

    complex(WP)         :: b(this%n_e*this%n)
    type(ext_complex_t) :: det
    complex(WP)         :: y_sh(this%n_e,this%n)
    logical             :: same_grid

    ! Reconstruct the solution on the shooting grid

    call this%build(omega)

    call this%sm%null_vector(b, det)

    y_sh = RESHAPE(b, SHAPE(y_sh))

    discrim = scale(det, -this%e_norm)    

    ! Build the recon grid

    this%recon_gp%omega_a = REAL(omega)
    this%recon_gp%omega_b = REAL(omega)

    call build_grid(this%recon_gp, this%bc, this%op, this%x, x)

    if(SIZE(x) == SIZE(this%x)) then
       same_grid = ALL(x == this%x)
    else
       same_grid = .FALSE.
    endif

    ! Reconstruct the full solution

    if(same_grid) then

       y = y_sh

    else

       allocate(y(this%n_e,SIZE(x)))

       call this%sh%recon(omega, this%x, y_sh, x, y, this%x_ad)

    endif

    ! Finish

    return

  end subroutine recon

!****

  function mode (this, omega) result (md)

    class(nad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega
    type(mode_t)                    :: md

    real(WP), allocatable    :: x(:)
    complex(WP), allocatable :: y(:,:)
    type(ext_complex_t)      :: discrim

    ! Reconstruct the solution

    call this%recon(omega, x, y, discrim)

    ! Initialize the mode
    
    call md%init(this%bc, this%op, omega, x, y, discrim, this%tc)

    ! Finish

    return

  end function mode

end module gyre_nad_bvp
