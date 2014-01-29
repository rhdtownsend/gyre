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
  use gyre_cocache
  $if ($MPI)
  use gyre_coeffs_mpi
  $endif
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_discfunc
  use gyre_shooter_nad
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

  type, extends (bvp_t) :: bvp_nad_t
     private
     class(coeffs_t), pointer       :: cf => null()
     type(cocache_t)                :: cc
     class(jacobian_t), allocatable :: jc
     class(ivp_t), allocatable      :: iv
     class(ivp_t), allocatable      :: iv_upw
     class(bound_t), allocatable    :: bd
     type(shooter_nad_t)            :: sh
     type(sysmtx_t)                 :: sm
     type(oscpar_t)                 :: op
     type(numpar_t)                 :: np
     type(gridpar_t), allocatable   :: shoot_gp(:)
     type(gridpar_t), allocatable   :: recon_gp(:)
     real(WP), allocatable          :: x_in(:)
     real(WP), allocatable          :: x(:)
     real(WP)                       :: x_upw
     integer, public                :: n
     integer, public                :: n_e
   contains 
     private
     $if ($GFORTRAN_PR57922)
     procedure, public :: final => final_
     $endif
     procedure, public :: set_x_upw => set_x_upw_
     procedure, public :: discrim => discrim_
     procedure         :: build_
     procedure         :: recon_
     procedure, public :: mode => mode_
     procedure, public :: coeffs => coeffs_
  end type bvp_nad_t

  ! Interfaces

  interface bvp_nad_t
     module procedure bvp_nad_t_
  end interface bvp_nad_t

  $if ($MPI)
  interface bcast
     module procedure bcast_
  end interface bcast
  $endif

  ! Access specifiers

  private

  public :: bvp_nad_t
  $if ($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  function bvp_nad_t_ (cf, op, np, shoot_gp, recon_gp, x_in) result (bp)

    use gyre_jacobian_nad_dziem
    use gyre_jacobian_nad_jcd

    use gyre_bound_nad_zero
    use gyre_bound_nad_dziem
    use gyre_bound_nad_unno
    use gyre_bound_nad_jcd

    use gyre_ivp_magnus_GL2
    use gyre_ivp_magnus_GL4
    use gyre_ivp_magnus_GL6
    use gyre_ivp_colloc_GL2
    use gyre_ivp_colloc_GL4
    use gyre_ivp_findiff_upw

    class(coeffs_t), pointer, intent(in) :: cf
    type(oscpar_t), intent(in)           :: op
    type(numpar_t), intent(in)           :: np
    type(gridpar_t), intent(in)          :: shoot_gp(:)
    type(gridpar_t), intent(in)          :: recon_gp(:)
    real(WP), allocatable, intent(in)    :: x_in(:)
    type(bvp_nad_t), target              :: bp

    integer               :: n
    real(WP), allocatable :: x_cc(:)

    ! Construct the bvp_nad_t

    ! Store parameters

    bp%op = op
    bp%np = np

    bp%shoot_gp = shoot_gp
    bp%recon_gp = recon_gp

    ! Set up the coefficient pointer
    
    bp%cf => cf

    ! Initialize the jacobian

    select case (bp%op%variables_type)
    case ('DZIEM')
       allocate(bp%jc, SOURCE=jacobian_nad_dziem_t(bp%cf, bp%op))
    case ('JCD')
       allocate(bp%jc, SOURCE=jacobian_nad_jcd_t(bp%cf, bp%op))
    case default
       $ABORT(Invalid variables_type)
    end select

    ! Initialize the boundary conditions

    select case (bp%op%outer_bound_type)
    case ('ZERO')
       allocate(bp%bd, SOURCE=bound_nad_zero_t(bp%cf, bp%jc, bp%op))
    case ('DZIEM')
       allocate(bp%bd, SOURCE=bound_nad_dziem_t(bp%cf, bp%jc, bp%op))
    case ('UNNO')
       allocate(bp%bd, SOURCE=bound_nad_unno_t(bp%cf, bp%jc, bp%op))
    case ('JCD')
       allocate(bp%bd, SOURCE=bound_nad_jcd_t(bp%cf, bp%jc, bp%op))
    case default
       $ABORT(Invalid bound_type)
    end select

    ! Initialize the IVP solvers

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

    allocate(bp%iv_upw, SOURCE=ivp_findiff_upw_t(bp%jc))

    select type (iv => bp%iv_upw)
    class is (ivp_findiff_upw_t)
       iv%stencil = [0,0,0,0,1,-1]
    end select

    ! Initialize the shooter

    bp%sh = shooter_nad_t(bp%cf, bp%iv, bp%iv_upw, bp%op, bp%np)

    ! Build the shooting grid

    call build_grid(bp%shoot_gp, bp%cf, bp%op, x_in, bp%x)

    n = SIZE(bp%x)

    ! Initialize the system matrix

    bp%sm = sysmtx_t(n-1, bp%sh%n_e, bp%bd%n_i, bp%bd%n_o)

    ! Other stuff

    if(ALLOCATED(x_in)) bp%x_in = x_in

    bp%x_upw = 0._WP

    bp%n = n
    bp%n_e = bp%sh%n_e

    ! Set up the coefficient cache

    x_cc = [bp%x(1),bp%sh%abscissa(bp%x),bp%x(n)]

    call bp%cf%attach_cache(bp%cc)
    call bp%cf%fill_cache(x_cc)
    call bp%cf%detach_cache()

    ! Finish

    return

  end function bvp_nad_t_

!****

  $if ($GFORTRAN_PR57922)

  subroutine final_ (this)

    class(bvp_nad_t), intent(inout) :: this

    ! Finalize the bvp_nad_t

    deallocate(this%jc)
    deallocate(this%iv)
    deallocate(this%bd)
    
    ! Finish

    return

  end subroutine final_

  $endif

!****

  subroutine set_x_upw_ (this, omega)

    class(bvp_nad_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega

    integer :: k

    ! Decide where to switch from the upwinded equations (interior)
    ! to the non-upwinded ones ones (exterior)

    this%x_upw = 0._WP

    x_upw_loop : do k = this%n,2,-1

       if(this%cf%tau_thm(this%x(k))*REAL(omega)*this%np%theta_ad > 1._WP) then
          this%x_upw = this%x(k)
          exit x_upw_loop
       endif

    end do x_upw_loop

    if(this%x_upw /= 0._WP) print *,'Set x_upw:',this%x_upw

!    this%x_upw = 0.

    ! Finish

    return

  end subroutine set_x_upw_

!****

  function discrim_ (this, omega, use_real)

    class(bvp_nad_t), intent(inout) :: this
    complex(WP), intent(in)         :: omega
    logical, optional, intent(in)   :: use_real
    type(ext_complex_t)             :: discrim_

    ! Evaluate the discriminant as the determinant of the sysmtx

    call this%build_(omega)

    call this%sm%determinant(discrim_, use_real, this%np%use_banded)

    ! Finish

    return

  end function discrim_

!****

  subroutine build_ (this, omega)

    class(bvp_nad_t), target, intent(inout) :: this
    complex(WP), intent(in)                 :: omega

    ! Set up the sysmtx

    call this%cf%attach_cache(this%cc)

    call this%sm%set_inner_bound(this%bd%inner_bound(this%x(1), omega), ext_complex_t(1._WP))
    call this%sm%set_outer_bound(this%bd%outer_bound(this%x(this%n), omega), ext_complex_t(1._WP))
 
    call this%sh%shoot(omega, this%x, this%sm, this%x_upw)

    call this%cf%detach_cache()

    call this%sm%scale_rows()

    ! Finish

    return

  end subroutine build_

!****

  subroutine recon_ (this, omega, x, y, x_ref, y_ref, discrim)

    class(bvp_nad_t), intent(inout)       :: this
    complex(WP), intent(in)               :: omega
    real(WP), allocatable, intent(out)    :: x(:)
    complex(WP), allocatable, intent(out) :: y(:,:)
    real(WP), intent(out)                 :: x_ref
    complex(WP), intent(out)              :: y_ref(:)
    type(ext_complex_t), intent(out)      :: discrim

    complex(WP) :: b(this%n_e*this%n)
    complex(WP) :: y_sh(this%n_e,this%n)
    logical     :: same_grid
    complex(WP) :: y_ref_(this%n_e,1)

    $CHECK_BOUNDS(SIZE(y_ref),this%n_e)

    ! Reconstruct the solution on the shooting grid

    call this%build_(omega)

    call this%sm%null_vector(b, discrim, this%np%use_banded)

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

       call this%sh%recon(omega, this%x, y_sh, x, y, this%x_upw)

    endif

    ! Reconstruct the solution at x_ref
    
    x_ref = MIN(MAX(this%op%x_ref, this%x(1)), this%x(this%n))
    
    call this%sh%recon(omega, this%x, y_sh, [x_ref], y_ref_)

    y_ref = y_ref_(:,1)

    ! Finish

    return

  end subroutine recon_

!****

  function mode_ (this, omega, discrim, use_real, omega_def) result (md)

    class(bvp_nad_t), target, intent(inout)   :: this
    complex(WP), intent(in)                   :: omega(:)
    type(ext_complex_t), optional, intent(in) :: discrim(:)
    logical, optional, intent(in)             :: use_real
    complex(WP), optional, intent(in)         :: omega_def(:)
    type(mode_t)                              :: md

    logical                   :: use_real_
    type(ext_complex_t)       :: omega_a
    type(ext_complex_t)       :: omega_b
    type(ext_complex_t)       :: discrim_a
    type(ext_complex_t)       :: discrim_b
    type(ext_real_t)          :: discrim_norm
    type(discfunc_t)          :: df
    integer                   :: n_iter_def
    integer                   :: n_iter
    complex(WP)               :: omega_root
    real(WP), allocatable     :: x(:)
    complex(WP), allocatable  :: y(:,:)
    real(WP)                  :: x_ref
    complex(WP)               :: y_ref(this%n_e)
    type(ext_complex_t)       :: discrim_root
    integer                   :: n
    integer                   :: i
    complex(WP), allocatable  :: y_c(:,:)
    complex(WP)               :: y_c_ref(6)
    type(ext_real_t)          :: chi 
    
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

    omega_a = ext_complex_t(omega(1))
    omega_b = ext_complex_t(omega(2))

    call this%set_x_upw(cmplx(0.5_WP*(omega_a+omega_b)))

    if(PRESENT(discrim)) then
       discrim_a = discrim(1)
       discrim_b = discrim(2)
    else
       discrim_a = this%discrim(cmplx(omega_a))
       discrim_b = this%discrim(cmplx(omega_b))
    endif
    
    discrim_norm = MAX(ABS(discrim_a), ABS(discrim_b))

    ! Set up the discriminant function

    df = discfunc_t(this)

    ! If omega_def is provided, do a preliminary root find using the
    ! deflated discriminant

    if(.FALSE. .AND. PRESENT(omega_def)) then

       ! (Don't pass discrim_a and discrim_b in/out, because they
       ! haven't been deflated)

       df%omega_def = omega_def

       n_iter_def = this%np%n_iter_max

       call df%narrow_pair(omega_a, omega_b, ext_real_t(0._WP), n_iter=n_iter_def)

       $ASSERT(n_iter_def <= this%np%n_iter_max,Too many iterations)

       deallocate(df%omega_def)

       ! If necessary, reset omega_a and omega_b so they are not
       ! coincident

       if(omega_b == omega_a) then
          omega_b = omega_a + TINY(0._WP)*(omega_a/ABS(omega_a))
       endif

       call df%expand_pair(omega_a, omega_b, ext_real_t(0._WP), discrim_a, discrim_b)

    endif

    ! Find the discriminant root

    n_iter = this%np%n_iter_max
 
   if(use_real_) then
       omega_root = real(df%root(ext_real_t(omega_a), ext_real_t(omega_b), ext_real_t(0._WP), &
                                 f_ex_a=ext_real_t(discrim_a), f_ex_b=ext_real_t(discrim_b), n_iter=n_iter))
    else
       omega_root = cmplx(df%root(omega_a, omega_b, ext_real_t(0._WP), &
                                  f_ez_a=discrim_a, f_ez_b=discrim_b, n_iter=n_iter))
    endif

    $ASSERT(n_iter <= this%np%n_iter_max,Too many iterations)

    ! Reconstruct the solution

    call this%recon_(omega_root, x, y, x_ref, y_ref, discrim_root)

    ! Calculate canonical variables

    n = SIZE(x)

    allocate(y_c(6,n))

    !$OMP PARALLEL DO 
    do i = 1,n
       y_c(:,i) = MATMUL(this%jc%trans_matrix(x(i), omega_root, .TRUE.), y(:,i))
    end do

    y_c_ref = MATMUL(this%jc%trans_matrix(x_ref, omega_root, .TRUE.), y_ref)

    ! Initialize the mode
    
    chi = ABS(discrim_root)/discrim_norm
    
    if(PRESENT(omega_def)) then
       md = mode_t(this%cf, this%op, omega_root, x, y_c, x_ref, y_c_ref, chi, n_iter)
    else
       md = mode_t(this%cf, this%op, omega_root, x, y_c, x_ref, y_c_ref, chi, n_iter)
    endif

    ! Finish

    return

  end function mode_

!****

  function coeffs_ (this) result (cf)

    class(bvp_nad_t), intent(in) :: this
    class(coeffs_t), pointer     :: cf

    ! Return the coefficients pointer

    cf => this%cf

    ! Finish

    return

  end function coeffs_

!****

  $if ($MPI)

  subroutine bcast_ (bp, root_rank, cf)

    type(bvp_nad_t), intent(inout)      :: bp
    integer, intent(in)                 :: root_rank
    class(coeffs_t), intent(in), target :: cf

    type(oscpar_t)               :: op
    type(numpar_t)               :: np
    type(gridpar_t), allocatable :: shoot_gp(:)
    type(gridpar_t), allocatable :: recon_gp(:)
    real(WP), allocatable        :: x_in(:)

    ! Broadcast the bvp_nad_t

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

       bp = bvp_nad_t(cf, op, np, shoot_gp, recon_gp, x_in)

    endif

    ! Finish

    return

  end subroutine bcast_

  $endif

end module gyre_bvp_nad
