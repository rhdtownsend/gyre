! Module   : gyre_bvp_rad
! Purpose  : solve radial adiabatic BVPs
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

module gyre_bvp_rad

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_bvp
  use gyre_coeffs
  use gyre_cocache
  $if($MPI)
  use gyre_coeffs_mpi
  $endif
  use gyre_oscpar
  use gyre_gridpar
  use gyre_numpar
  use gyre_discfunc
  use gyre_shooter_rad
  use gyre_jacobian
  use gyre_ivp
  use gyre_bound
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_grid
  use gyre_mode

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(bvp_t) :: bvp_rad_t
     private
     class(coeffs_t), pointer       :: cf => null()
     type(cocache_t)                :: cc
     class(jacobian_t), allocatable :: jc
     class(ivp_t), allocatable      :: iv
     class(bound_t), allocatable    :: bd
     type(shooter_rad_t)            :: sh
     type(sysmtx_t)                 :: sm
     type(oscpar_t)                 :: op
     type(numpar_t)                 :: np
     type(gridpar_t), allocatable   :: shoot_gp(:)
     type(gridpar_t), allocatable   :: recon_gp(:)
     real(WP), allocatable          :: x_in(:)
     real(WP), allocatable          :: x(:)
     integer, public                :: n
     integer, public                :: n_e
   contains 
     private
     $if($GFORTRAN_PR57922)
     procedure, public :: final
     $endif
     procedure, public :: discrim
     procedure         :: build
     procedure         :: recon
     procedure, public :: mode
     procedure, public :: coeffs
  end type bvp_rad_t

  ! Interfaces

  interface bvp_rad_t
     module procedure init_bp
  end interface bvp_rad_t

  $if ($MPI)
  interface bcast
     module procedure bcast_bp
  end interface bcast
  $endif

  ! Access specifiers

  private

  public :: bvp_rad_t
  $if ($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  function init_bp (cf, op, np, shoot_gp, recon_gp, x_in) result (bp)

    use gyre_jacobian_rad_dziem
    use gyre_jacobian_rad_jcd
    use gyre_jacobian_rad_mix

    use gyre_bound_rad_zero
    use gyre_bound_rad_dziem
    use gyre_bound_rad_unno
    use gyre_bound_rad_jcd

    use gyre_ivp_magnus_GL2
    use gyre_ivp_magnus_GL4
    use gyre_ivp_magnus_GL6
    use gyre_ivp_colloc_GL2
    use gyre_ivp_colloc_GL4

    class(coeffs_t), pointer, intent(in) :: cf
    type(oscpar_t), intent(in)           :: op
    type(numpar_t), intent(in)           :: np
    type(gridpar_t), intent(in)          :: shoot_gp(:)
    type(gridpar_t), intent(in)          :: recon_gp(:)
    real(WP), allocatable, intent(in)    :: x_in(:)
    type(bvp_rad_t), target              :: bp

    integer               :: n
    real(WP), allocatable :: x_cc(:)

    $ASSERT(op%l == 0,Invalid harmonic degree)

    ! Construct the bvp_rad

    ! Store parameters

    bp%op = op
    bp%np = np

    bp%shoot_gp = shoot_gp
    bp%recon_gp = recon_gp

    ! Set up the coefficient pointer
    
    bp%cf => cf

    ! Initialize the jacobian

    select case (op%variables_type)
    case ('DZIEM')
       allocate(bp%jc, SOURCE=jacobian_rad_dziem_t(bp%cf, bp%op))
    case ('JCD')
       allocate(bp%jc, SOURCE=jacobian_rad_jcd_t(bp%cf, bp%op))
    case ('MIX')
       allocate(bp%jc, SOURCE=jacobian_rad_mix_t(bp%cf, bp%op))
    case default
       $ABORT(Invalid variables_type)
    end select

    ! Initialize the boundary conditions

    select case (bp%op%outer_bound_type)
    case ('ZERO')
       allocate(bp%bd, SOURCE=bound_rad_zero_t(bp%cf, bp%jc, bp%op))
    case ('DZIEM')
       allocate(bp%bd, SOURCE=bound_rad_dziem_t(bp%cf, bp%jc, bp%op))
    case ('UNNO')
       allocate(bp%bd, SOURCE=bound_rad_unno_t(bp%cf, bp%jc, bp%op))
    case ('JCD')
       allocate(bp%bd, SOURCE=bound_rad_jcd_t(bp%cf, bp%jc, bp%op))
    case default
       $ABORT(Invalid bound_type)
    end select

    ! Initialize the IVP solver

    select case (bp%np%ivp_solver_type)
    case ('MAGNUS_GL2')
       allocate(bp%iv, SOURCE=ivp_magnus_GL2_t(bp%jc))
    case ('MAGNUS_GL4')
       allocate(bp%iv, SOURCE=ivp_magnus_GL4_t(bp%jc))
    case ('MAGNUS_GL6')
       allocate(bp%iv, SOURCE=ivp_magnus_GL6_t(bp%jc))
    case ('FINDIFF_GL2')
       allocate(bp%iv, SOURCE=ivp_colloc_GL2_t(bp%jc))
    case ('FINDIFF_GL4')
       allocate(bp%iv, SOURCE=ivp_colloc_GL4_t(bp%jc))
    case default
       $ABORT(Invalid ivp_solver_type)
    end select

    ! Initialize the shooter

    bp%sh = shooter_rad_t(bp%cf, bp%iv, bp%op, bp%np)

    ! Build the shooting grid

    call build_grid(bp%shoot_gp, bp%cf, bp%op, x_in, bp%x, verbose=.TRUE.)

    n = SIZE(bp%x)

    ! Initialize the system matrix

    bp%sm = sysmtx_t(n-1, bp%jc%n_e, bp%bd%n_i, bp%bd%n_o)

    ! Other stuff

    if(ALLOCATED(x_in)) bp%x_in = x_in

    bp%n = n
    bp%n_e = bp%sh%n_e

    ! Set up the coefficient cache

    x_cc = [bp%x(1),bp%sh%abscissa(bp%x),bp%x(n)]

    call bp%cf%attach_cache(bp%cc)
    call bp%cf%fill_cache(x_cc)
    call bp%cf%detach_cache()

    ! Finish

    return

  end function init_bp

!****

  $if ($GFORTRAN_PR57922)

  subroutine final (this)

    class(bvp_rad_t), intent(inout) :: this

    ! Finalize the bvp_rad

    call this%cc%final()
    call this%sm%final()

    deallocate(this%jc)
    deallocate(this%iv)
    deallocate(this%bd)

    deallocate(this%shoot_gp)
    deallocate(this%recon_gp)

    deallocate(this%x)
    deallocate(this%x_in)
    
    ! Finish

    return

  end subroutine final

  $endif

!****

  $if ($MPI)

  subroutine bcast_bp (this, root_rank, cf)

    class(bvp_rad_t), intent(inout)     :: bp
    integer, intent(in)                 :: root_rank
    class(coeffs_t), intent(in), target :: cf

    type(oscpar_t)                :: op
    type(numpar_t)                :: np
    type(gridpar_t), allocatable  :: shoot_gp(:)
    type(gridpar_t), allocatable  :: recon_gp(:)
    real(WP), allocatable         :: x_in(:)

    ! Broadcast the bvp

    if(MPI_RANK == root_rank) then

       call bcast(bp%op, root_rank)
       call bcast(bp%np, root_rank)

       call bcast_alloc(bp%shoot_gp, root_rank)
       call bcast_alloc(bp%recon_gp, root_rank)

       call bcast_alloc(bp%x_in, root_rank)

    else

       call bcast(op, root_rank)
       call bcast(np, root_rank)

       call bcast_alloc(shoot_gp, root_rank)
       call bcast_alloc(recon_gp, root_rank)

       call bcast_alloc(x_in, root_rank)

       call bp%init(cf, op, np, shoot_gp, recon_gp, x_in)

    endif

    ! Finish

    return

  end subroutine bcast_bp

  $endif

!****

  function discrim (this, omega, use_real)

    class(bvp_rad_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega
    logical, intent(in), optional   :: use_real
    type(ext_complex_t)             :: discrim

    ! Evaluate the discriminant as the determinant of the sysmtx

    call this%build(omega)

    call this%sm%determinant(discrim, use_real, this%np%use_banded)

    ! Finish

    return

  end function discrim

!****

  subroutine build (this, omega)

    class(bvp_rad_t), target, intent(inout) :: this
    complex(WP), intent(in)                 :: omega

    ! Set up the sysmtx

    call this%cf%attach_cache(this%cc)

    call this%sm%set_inner_bound(this%bd%inner_bound(this%x(1), omega), ext_complex(1._WP))
    call this%sm%set_outer_bound(this%bd%outer_bound(this%x(this%n), omega), ext_complex(1._WP))

    call this%sh%shoot(omega, this%x, this%sm)

    call this%cf%detach_cache()

    call this%sm%scale_rows()

    ! Finish

    return

  end subroutine build

!****

  subroutine recon (this, omega, x, y, x_ref, y_ref, discrim, use_real)

    class(bvp_rad_t), intent(inout)       :: this
    complex(WP), intent(in)               :: omega
    real(WP), allocatable, intent(out)    :: x(:)
    complex(WP), allocatable, intent(out) :: y(:,:)
    real(WP), intent(out)                 :: x_ref
    complex(WP), intent(out)              :: y_ref(:)
    type(ext_complex_t), intent(out)      :: discrim
    logical, intent(in), optional         :: use_real

    complex(WP) :: b(this%n_e*this%n)
    complex(WP) :: y_sh(this%n_e,this%n)
    logical     :: same_grid
    complex(WP) :: y_ref_(this%n_e,1)

    $CHECK_BOUNDS(SIZE(y_ref),this%n_e)

    ! Reconstruct the solution on the shooting grid

    call this%build(omega)

    call this%sm%null_vector(b, discrim, use_real, this%np%use_banded)

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

       call this%sh%recon(omega, this%x, y_sh, x, y)

    endif

    ! Reconstruct the solution at x_ref
    
    x_ref = MIN(MAX(this%op%x_ref, this%x(1)), this%x(this%n))

    call this%sh%recon(omega, this%x, y_sh, [x_ref], y_ref_)

    y_ref = y_ref_(:,1)

    ! Finish

    return

  end subroutine recon

!****

  function mode (this, omega, discrim, use_real, omega_def) result (md)

    class(bvp_rad_t), target, intent(inout)   :: this
    complex(WP), intent(in)                   :: omega(:)
    type(ext_complex_t), intent(in), optional :: discrim(:)
    logical, intent(in), optional             :: use_real
    complex(WP), intent(in), optional         :: omega_def(:)
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
    real(WP)                 :: x_ref
    complex(WP)              :: y_ref(this%n_e)
    type(ext_complex_t)      :: discrim_root
    integer                  :: n
    integer                  :: i
    complex(WP), allocatable :: y_c(:,:)
    complex(WP)              :: y_c_ref(6)
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

    ! Set up the discriminant function

    df = discfunc_t(this)

    if(PRESENT(omega_def)) df%omega_def = omega_def

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

    call this%recon(omega_root, x, y, x_ref, y_ref, discrim_root)

    ! Calculate canonical variables

    n = SIZE(x)

    allocate(y_c(6,n))

    !$OMP PARALLEL DO 
    do i = 1,n
       y_c(1:2,i) = MATMUL(this%jc%trans_matrix(x(i), omega_root, .TRUE.), y(:,i))
       y_c(3,i) = 0._WP
       y_c(4,i) = -y_c(1,i)*this%cf%U(x(i))
       y_c(5:6,i) = 0._WP
    end do

    y_c_ref(1:2) = MATMUL(this%jc%trans_matrix(x_ref, omega_root, .TRUE.), y_ref)
    y_c_ref(3) = 0._WP
    y_c_ref(4) = -y_c_ref(1)*this%cf%U(x_ref)
    y_c_ref(5:6) = 0._WP

    ! Initialize the mode
    
    chi = ABS(discrim_root)/MAX(ABS(discrim_a), ABS(discrim_b))
    
    md = mode_t(this%cf, this%op, omega_root, x, y_c, x_ref, y_c_ref, chi, n_iter)

    ! Finish

    return

  end function mode

!****

  function coeffs (this) result (cf)

    class(bvp_rad_t), intent(in) :: this
    class(coeffs_t), pointer     :: cf

    ! Return the coefficients pointer

    cf => this%cf

    ! Finish

    return

  end function coeffs

end module gyre_bvp_rad
