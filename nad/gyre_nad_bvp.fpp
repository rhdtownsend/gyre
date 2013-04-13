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
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_grid
  use gyre_eigfunc

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
     procedure, public :: discrim
     procedure         :: build
     procedure         :: recon
     procedure, public :: eigfunc
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

  subroutine init (this, bc, tc, op, np, shoot_gp, recon_gp, x_in)

    class(nad_bvp_t), intent(out)     :: this
    class(base_coeffs_t), intent(in)  :: bc
    class(therm_coeffs_t), intent(in) :: tc
    type(oscpar_t), intent(in)        :: op
    type(numpar_t), intent(in)        :: np
    type(gridpar_t), intent(in)       :: shoot_gp(:)
    type(gridpar_t), intent(in)       :: recon_gp(:)
    real(WP), allocatable, intent(in) :: x_in(:)

    integer :: n

    ! Initialize the nad_bvp

    ! Create the shooting grid

    call build_grid(shoot_gp, bc, op, x_in, this%x)

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

    this%e_norm = 0

    this%n = n
    this%n_e = this%sh%n_e

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

    discrim = this%sm%determinant()

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

    call this%sm%set_inner_bound(this%bd%inner_bound(omega))
    call this%sm%set_outer_bound(this%bd%outer_bound(omega))

    call this%sh%shoot(omega, this%x, this%sm)

    ! Finish

    return

  end subroutine build

!****

  subroutine recon (this, omega, x, y)

    class(nad_bvp_t), intent(inout)       :: this
    complex(WP), intent(in)               :: omega
    real(WP), allocatable, intent(out)    :: x(:)
    complex(WP), allocatable, intent(out) :: y(:,:)

    complex(WP) :: y_sh(this%n_e,this%n)

    ! Reconstruct the solution on the shooting grid

    call this%build(omega)

    y_sh = RESHAPE(this%sm%null_vector(), SHAPE(y_sh))

    ! Reconstruct the solution on the shooting grid

    call this%build(omega)

    y_sh = RESHAPE(this%sm%null_vector(), SHAPE(y_sh))

    ! Build the recon grid

    this%recon_gp%omega_a = REAL(omega)
    this%recon_gp%omega_b = REAL(omega)

    call build_grid(this%recon_gp, this%bc, this%op, this%x, x)

    ! Reconstruct the full solution

    allocate(y(this%n_e,SIZE(x)))

    call this%sh%recon(omega, this%x, y_sh, x, y)

    ! Finish

    return

    ! Finish

    return

  end subroutine recon

!****

  function eigfunc (this, omega) result (ef)

    class(nad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega
    type(eigfunc_t)                 :: ef

    real(WP), allocatable    :: x(:)
    complex(WP), allocatable :: y(:,:)

    ! Reconstruct the solution

    call this%recon(omega, x, y)

    ! Initialize the eigfunc
    
    call ef%init(this%bc, this%tc, this%op, omega, x, y)

    ! Finish

    return

  end function eigfunc

end module gyre_nad_bvp
