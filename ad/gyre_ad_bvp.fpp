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
  use core_parallel

  use gyre_bvp
  use gyre_base_coeffs
  use gyre_therm_coeffs
  $if($MPI)
  use gyre_base_coeffs_mpi
  use gyre_therm_coeffs_mpi
  $endif
  use gyre_oscpar
  use gyre_gridpar
  use gyre_numpar
  use gyre_ad_shooter
  use gyre_ad_bound
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_grid
  use gyre_mode

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(bvp_t) :: ad_bvp_t
     private
     class(base_coeffs_t), allocatable  :: bc
     class(therm_coeffs_t), allocatable :: tc
     type(oscpar_t)                     :: op
     type(numpar_t)                     :: np
     type(gridpar_t), allocatable       :: shoot_gp(:)
     type(gridpar_t), allocatable       :: recon_gp(:)
     type(ad_shooter_t)                 :: sh
     type(ad_bound_t)                   :: bd
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
!     procedure, public :: get_op
     procedure, public :: get_np
!     procedure, public :: get_shoot_gp
!     procedure, public :: get_recon_gp
     procedure, public :: set_norm
     procedure, public :: discrim
     procedure         :: build
     procedure         :: recon
     procedure, public :: mode
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

  subroutine init (this, bc, op, np, shoot_gp, recon_gp, x_in, tc)

    class(ad_bvp_t), intent(out)                :: this
    class(base_coeffs_t), intent(in)            :: bc
    type(oscpar_t), intent(in)                  :: op
    type(numpar_t), intent(in)                  :: np
    type(gridpar_t), intent(in)                 :: shoot_gp(:)
    type(gridpar_t), intent(in)                 :: recon_gp(:)
    real(WP), allocatable, intent(in)           :: x_in(:)
    class(therm_coeffs_t), intent(in), optional :: tc

    integer               :: n
    real(WP), allocatable :: x_cc(:)

    ! Initialize the ad_bvp

    ! Create the shooting grid

    call build_grid(shoot_gp, bc, op, x_in, this%x, tc)

    n = SIZE(this%x)

    ! Set up components
    
    allocate(this%bc, SOURCE=bc)
    if(PRESENT(tc)) allocate(this%tc, SOURCE=tc)

    this%op = op
    this%np = np

    this%shoot_gp = shoot_gp
    this%recon_gp = recon_gp

    call this%sh%init(this%bc, this%op, this%np)
    call this%bd%init(this%bc, this%op)

    call this%sm%init(n-1, this%sh%n_e, this%bd%n_i, this%bd%n_o)

    if(ALLOCATED(x_in)) this%x_in = x_in

    this%e_norm = 0

    this%n = n
    this%n_e = this%sh%n_e

    ! Set up the coefficient caches

    x_cc = [this%x(1),this%sh%abscissa(this%x),this%x(n)]

    call this%bc%fill_cache(x_cc)
    if(PRESENT(tc)) call this%tc%fill_cache(x_cc)

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_bp (this, root_rank)

    class(ad_bvp_t), intent(inout) :: this
    integer, intent(in)            :: root_rank

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

    class(ad_bvp_t), intent(in), target :: this
    $ITEM_TYPE, pointer                 :: item

    ! Get a pointer to the item

    item => this%$ITEM_NAME

    ! Finish

    return

  end function get_$ITEM_NAME
  
  $endsub

!  $GET(bc,class(base_coeffs_t))
!  $GET(op,type(oscpar_t))
  $GET(np,type(numpar_t))
!  $GET(shoot_gp,type(gridpar_t))
!  $GET(recon_gp,type(gridpar_t))
    
!****

  subroutine set_norm (this, omega)

    class(ad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)        :: omega

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

    discrim = this%sm%determinant(use_real=.TRUE.)

    ! Apply the normalization

    if(norm_) discrim = scale(discrim, -this%e_norm)

    ! Finish

    return

  end function discrim

!****

  subroutine build (this, omega)

    class(ad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)        :: omega

    ! Set up the sysmtx

    call this%bc%enable_cache()
    if(ALLOCATED(this%tc)) call this%tc%enable_cache()

    call this%sm%set_inner_bound(this%bd%inner_bound(this%x(1), omega))
    call this%sm%set_outer_bound(this%bd%outer_bound(this%x(this%n), omega))

    call this%sh%shoot(omega, this%x, this%sm)

    call this%bc%disable_cache()
    if(ALLOCATED(this%tc)) call this%tc%disable_cache()

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

  end subroutine recon

!****

  function mode (this, omega) result (md)

    class(ad_bvp_t), intent(inout) :: this
    complex(WP), intent(in)        :: omega
    type(mode_t)                   :: md

    real(WP), allocatable    :: x(:)
    complex(WP), allocatable :: y(:,:)
    integer                  :: n
    complex(WP), allocatable :: y_6(:,:)

    ! Reconstruct the solution

    call this%recon(omega, x, y)

    n = SIZE(x)

    allocate(y_6(6,n))

    y_6(1:4,:) = y

    if(ALLOCATED(this%tc)) then
       call recon_y_6(this, omega, x, y_6)
       call recon_y_5(this, omega, x, y_6)
    else
       y_6(5,:) = 0._WP
       y_6(6,:) = 0._WP
    end if

    ! Initialize the mode
    
    call md%init(this%bc, this%tc, this%op, omega, x, y_6)

    ! Finish

    return

  end function mode

!****

  subroutine recon_y_5 (bp, omega, x, y)

    class(ad_bvp_t), intent(in) :: bp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x(:)
    complex(WP), intent(inout)  :: y(:,:)

    complex(WP) :: dy_6_dx(SIZE(x))
    complex(WP) :: A_6(SIZE(x),6)

    $CHECK_BOUNDS(SIZE(y, 1),6)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Reconstruct y_5 (the entropy perturbation) using the
    ! energy equation

    dy_6_dx = deriv(x, y(6,:))
    
    associate(V => bp%bc%V(x), c_1 => bp%bc%c_1(x), &
              nabla_ad => bp%bc%nabla_ad(x), nabla => bp%tc%nabla(x), &
              c_rad => bp%tc%c_rad(x), dc_rad => bp%tc%dc_rad(x), &
              c_thm => bp%tc%c_thm(x), &
              c_eps_ad => bp%tc%c_eps_ad(x), c_eps_S => bp%tc%c_eps_S(x), &
              l => bp%op%l)

      ! Calculate Jacobian coefficients. (cf. gyre_nad_jacobian; A_6
      ! has been multiplied by V)

      A_6(:,1) = V*(l*(l+1)*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad)
      A_6(:,2) = V*(V*c_eps_ad - l*(l+1)*c_rad*(nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega**2)))
      A_6(:,3) = V*(l*(l+1)*nabla_ad/nabla*c_rad - V*c_eps_ad)
      A_6(:,4) = 0._WP
      A_6(:,5) = V*c_eps_S - l*(l+1)*c_rad/nabla - V*(0._WP,1._WP)*omega*c_thm
      A_6(:,6) = V*(-1._WP - l)

      ! Set up y_5

      where(x /= 0._WP)
         y(5,:) = (x*V*dy_6_dx - (A_6(:,1)*y(1,:) + A_6(:,2)*y(2,:) + &
                                  A_6(:,3)*y(3,:) + A_6(:,4)*y(4,:) + A_6(:,6)*y(6,:)))/A_6(:,5)
      elsewhere
         y(5,:) = 0._WP
      end where

    end associate

    ! Finish

    return

  contains

    function deriv (x, y) result (dy_dx)

      real(WP), intent(in)    :: x(:)
      complex(WP), intent(in) :: y(:)
      complex(WP)             :: dy_dx(SIZE(x))
      
      integer :: n
      integer :: i

      $CHECK_BOUNDS(SIZE(y),SIZE(x))

      ! Differentiate y(x) using centered finite differences

      n = SIZE(x)

      dy_dx(1) = (y(2) - y(1))/(x(2) - x(1))

      do i = 2,n-1
         dy_dx(i) = 0.5_WP*((y(i) - y(i-1))/(x(i) - x(i-1)) + &
                            (y(i+1) - y(i))/(x(i+1) - x(i)))
      end do

      dy_dx(n) = (y(n) - y(n-1))/(x(n) - x(n-1))

      ! Finish

      return

    end function deriv

  end subroutine recon_y_5

!****

  subroutine recon_y_6 (bp, omega, x, y)

    class(ad_bvp_t), intent(in) :: bp
    complex(WP), intent(in)     :: omega
    real(WP), intent(in)        :: x(:)
    complex(WP), intent(inout)  :: y(:,:)

    complex(WP) :: A_5(SIZE(x),6)

    $CHECK_BOUNDS(SIZE(y, 1),6)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Reconstruct y_6 (the luminosity perturbation) using the
    ! quasi-adiabatic diffusion equation

    associate(U => bp%bc%U(x), c_1 => bp%bc%c_1(x), &
              nabla_ad => bp%bc%nabla_ad(x), &
              c_rad => bp%tc%c_rad(x), &
              c_dif => bp%tc%c_dif(x), nabla => bp%tc%nabla(x), &
              l => bp%op%l)

      ! Calculate Jacobian coefficients. (cf. gyre_nad_jacobian; A_5
      ! has been implicitly divided by V)

      A_5(:,1) = nabla_ad*(U - c_1*omega**2) - 4._WP*(nabla_ad - nabla) + c_dif
      A_5(:,2) = l*(l+1)/(c_1*omega**2)*(nabla_ad - nabla) - c_dif
      A_5(:,3) = c_dif
      A_5(:,4) = nabla_ad
      A_5(:,5) = 0._WP
      A_5(:,6) = -nabla/c_rad

      ! Set up y_6

      y(6,:) = -(A_5(:,1)*y(1,:) + A_5(:,2)*y(2,:) + A_5(:,3)*y(3,:) + A_5(:,4)*y(4,:))/A_5(:,6)

    end associate

    ! Finish

    return

  end subroutine recon_y_6

end module gyre_ad_bvp
