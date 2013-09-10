! Module   : gyre_bvp_nad
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

module gyre_bvp_nad

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_bvp
  use gyre_coeffs
  $if($MPI)
  use gyre_coeffs_mpi
  $endif
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_discfunc
  use gyre_shooter_nad
  use gyre_jacobian
  use gyre_bound
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_grid
  use gyre_mode

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (bvp_t) :: bvp_nad_t
     private
     class(coeffs_t), allocatable   :: cf
     class(jacobian_t), allocatable :: jc
     class(bound_t), allocatable    :: bd
     type(shooter_nad_t)            :: sh
     type(sysmtx_t)                 :: sm
     type(oscpar_t)                 :: op
     type(numpar_t)                 :: np
     type(gridpar_t), allocatable   :: shoot_gp(:)
     type(gridpar_t), allocatable   :: recon_gp(:)
     real(WP), allocatable          :: x_in(:)
     real(WP), allocatable          :: x(:)
     real(WP)                       :: x_ad
     integer                        :: e_norm
     integer, public                :: n
     integer, public                :: n_e
   contains 
     private
     procedure, public :: init
     procedure, public :: set_x_ad
     procedure, public :: discrim
     procedure         :: build
     procedure         :: recon
     procedure, public :: mode
  end type bvp_nad_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_bp
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: bvp_nad_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init (this, cf, op, np, shoot_gp, recon_gp, x_in)

    use gyre_bvp_nad_initmods

    class(bvp_nad_t), intent(out)     :: this
    class(coeffs_t), intent(in)       :: cf
    type(oscpar_t), intent(in)        :: op
    type(numpar_t), intent(in)        :: np
    type(gridpar_t), intent(in)       :: shoot_gp(:)
    type(gridpar_t), intent(in)       :: recon_gp(:)
    real(WP), allocatable, intent(in) :: x_in(:)

    integer               :: n
    real(WP), allocatable :: x_cc(:)

    ! Initialize the bvp_nad

    ! Store parameters

    this%op = op
    this%np = np

    this%shoot_gp = shoot_gp
    this%recon_gp = recon_gp

    ! Copy coefficients

    allocate(this%cf, SOURCE=cf)

    ! Initialize the jacobian

    select case (this%op%variables_type)
    case ('DZIEM')
       allocate(jacobian_nad_dziem_t::this%jc)
    case ('JCD')
       allocate(jacobian_nad_jcd_t::this%jc)
    case default
       $ABORT(Invalid variables_type)
    end select

    call this%jc%init(this%cf, this%op)

    ! Initialize the boundary conditions

    select case (this%op%outer_bound_type)
    case ('ZERO')
       allocate(bound_nad_zero_t::this%bd)
    case ('DZIEM')
       allocate(bound_nad_dziem_t::this%bd)
    ! case ('UNNO')
    !    allocate(bound_nad_unno_t::this%bd)
    ! case ('JCD')
    !    allocate(bound_nad_jcd_t::this%bd)
    case default
       $ABORT(Invalid bound_type)
    end select

    call this%bd%init(this%cf, this%jc, this%op)

    ! Initialize the shooter

    call this%sh%init(this%cf, this%jc, this%op, this%np)

    ! Build the shooting grid

    call build_grid(this%shoot_gp, this%cf, this%op, x_in, this%x)

    n = SIZE(this%x)

    ! Initialize the system matrix

    call this%sm%init(n-1, this%sh%n_e, this%bd%n_i, this%bd%n_o)

    ! Other stuff

    if(ALLOCATED(x_in)) this%x_in = x_in

    this%x_ad = 0._WP

    this%e_norm = 0

    this%n = n
    this%n_e = this%sh%n_e

    ! Set up the coefficient caches

    x_cc = [this%x(1),this%sh%abscissa(this%x),this%x(n)]

    call this%cf%fill_cache(x_cc)

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_bp (this, root_rank)

    class(bvp_nad_t), intent(inout) :: this
    integer, intent(in)             :: root_rank

    class(coeffs_t), allocatable :: cf
    type(oscpar_t)               :: op
    type(numpar_t)               :: np
    type(gridpar_t), allocatable :: shoot_gp(:)
    type(gridpar_t), allocatable :: recon_gp(:)
    real(WP), allocatable        :: x_in(:)

    ! Broadcast the bvp

    if(MPI_RANK == root_rank) then

       call bcast_alloc(this%cf, root_rank)
      
       call bcast(this%op, root_rank)
       call bcast(this%np, root_rank)

       call bcast_alloc(this%shoot_gp, root_rank)
       call bcast_alloc(this%recon_gp, root_rank)

       call bcast_alloc(this%x_in, root_rank)

    else

       call bcast_alloc(cf, root_rank)

       call bcast(op, root_rank)
       call bcast(np, root_rank)

       call bcast_alloc(shoot_gp, root_rank)
       call bcast_alloc(recon_gp, root_rank)

       call bcast_alloc(x_in, root_rank)

       call this%init(cf, op, np, shoot_gp, recon_gp, x_in)

    endif

    ! Finish

    return

  end subroutine bcast_bp

  $endif

!****

  subroutine set_x_ad (this, omega)

    class(bvp_nad_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega

    integer :: k

    ! Decide where to switch from the adiabatic equations (interior)
    ! to the non-adiabastic ones (exterior)

    this%x_ad = 0._WP

    x_ad_loop : do k = this%n,2,-1

       if(this%cf%tau_thm(this%x(k))*REAL(omega)*this%np%theta_ad > 1._WP) then
          this%x_ad = this%x(k)
          exit x_ad_loop
       endif

    end do x_ad_loop

    ! Finish

    return

  end subroutine set_x_ad

!****

  function discrim (this, omega, use_real)

    class(bvp_nad_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega
    logical, intent(in), optional   :: use_real
    type(ext_complex_t)             :: discrim

    ! Evaluate the discriminant as the determinant of the sysmtx

    call this%build(omega)

    call this%sm%determinant(discrim, use_real, this%np%use_banded)

    ! Scale the discriminant using the normalizing exponent

    discrim = scale(discrim, -this%e_norm)

    ! Finish

    return

  end function discrim

!****

  subroutine build (this, omega)

    class(bvp_nad_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega

    ! Set up the sysmtx

    call this%cf%enable_cache()

    call this%sm%set_inner_bound(this%bd%inner_bound(this%x(1), omega))
    call this%sm%set_outer_bound(this%bd%outer_bound(this%x(this%n), omega))

    call this%sh%shoot(omega, this%x, this%sm, this%x_ad)

    call this%cf%disable_cache()

    ! Finish

    return

  end subroutine build

!****

  subroutine recon (this, omega, x, y, discrim)

    class(bvp_nad_t), intent(inout)       :: this
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

    call this%sm%null_vector(b, det, this%np%use_banded)

    discrim = scale(det, -this%e_norm)    

    y_sh = RESHAPE(b, SHAPE(y_sh))

    ! Build the recon grid

    this%recon_gp%omega_a = REAL(omega)
    this%recon_gp%omega_b = REAL(omega)

    call build_grid(this%recon_gp, this%cf, this%op, this%x, x)

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

  function mode (this, omega, discrim, use_real) result (md)

    class(bvp_nad_t), intent(inout)           :: this
    complex(WP), intent(in)                   :: omega(:)
    type(ext_complex_t), intent(in), optional :: discrim(:)
    logical, intent(in), optional             :: use_real
    type(mode_t)                              :: md

    logical                  :: use_real_
    complex(WP)              :: omega_a
    complex(WP)              :: omega_b
    type(ext_complex_t)      :: discrim_a
    type(ext_complex_t)      :: discrim_b
    type(discfunc_t)         :: df
    integer                  :: n_iter
    complex(WP)              :: omega_root
    real(WP), allocatable    :: x(:)
    complex(WP), allocatable :: y(:,:)
    type(ext_complex_t)      :: discrim_root
    integer                  :: n
    integer                  :: i
    complex(WP), allocatable :: y_c(:,:)
    type(ext_real_t)         :: chi 
    
    $CHECK_BOUNDS(SIZE(omega),2)
    
    if(PRESENT(discrim)) then
       $CHECK_BOUNDS(SIZE(discrim),2)
    endif

    if(PRESENT(use_real)) then
       use_real_ = use_real
    else
       use_real_ = .FALSE.
    endif

    ! Unpack arguments

    omega_a = omega(1)
    omega_b = omega(2)

    if(PRESENT(discrim)) then
       discrim_a = discrim(1)
       discrim_b = discrim(2)
    else
       discrim_a = this%discrim(omega_a)
       discrim_b = this%discrim(omega_b)
    endif

    ! Set the normalizing exponent

    this%e_norm = MAX(exponent(discrim_a), exponent(discrim_b))

    discrim_a = scale(discrim_a, -this%e_norm)
    discrim_b = scale(discrim_b, -this%e_norm)

    ! Set up the discriminant function

    call df%init(this)

    ! Find the discriminant root

    n_iter = this%np%n_iter_max

    if(use_real_) then
       omega_root = real(df%root(ext_real(omega_a), ext_real(omega_b), ext_real(0._WP), &
                            f_ex_a=ext_real(discrim_a), f_ex_b=ext_real(discrim_b), n_iter=n_iter))
    else
       omega_root = cmplx(df%root(ext_complex(omega_a), ext_complex(omega_b), ext_real(0._WP), &
                            f_ez_a=discrim_a, f_ez_b=discrim_b, n_iter=n_iter))
    endif

    $ASSERT(n_iter <= this%np%n_iter_max,Too many iterations)

    ! Reconstruct the solution

    call this%recon(omega_root, x, y, discrim_root)

    ! Calculate canonical variables

    n = SIZE(x)

    allocate(y_c(6,n))

    !$OMP PARALLEL DO 
    do i = 1,n
       y_c(:,i) = MATMUL(this%jc%trans_matrix(x(i), omega_root, .TRUE.), y(:,i))
    end do

    ! Initialize the mode
    
    chi = ABS(discrim_root)/MAX(ABS(discrim_a), ABS(discrim_b))
    
    call md%init(this%cf, this%op, omega_root, x, y_c, chi, n_iter)

    ! Reset the normalizing exponent

    this%e_norm = 0

    ! Finish

    return

  end function mode

end module gyre_bvp_nad
