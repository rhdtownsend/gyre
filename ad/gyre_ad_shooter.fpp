! Module   : gyre_ad_shooter
! Purpose  : adiabatic multiple shooting
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

module gyre_ad_shooter

  ! Uses

  use core_kinds

  use gyre_mech_coeffs
  use gyre_oscpar
  use gyre_ad_jacobian
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_ivp
  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ad_shooter_t
     private
     class(mech_coeffs_t), pointer :: mc => null()
     type(oscpar_t), pointer       :: op => null()
     type(ad_jacobian_t)           :: jc
     real(WP), allocatable         :: x(:)
     real(WP)                      :: alpha_osc
     real(WP)                      :: alpha_exp
     integer                       :: n_center
     integer                       :: n_floor
     integer, public               :: n
     integer, public               :: n_e
     character(LEN=256)            :: solver_type
   contains
     private
     procedure, public :: init
     procedure, public :: shoot
     procedure, public :: recon => recon_sh
  end type ad_shooter_t

  ! Access specifiers

  private

  public :: ad_shooter_t

  ! Procedures

contains

  subroutine init (this, mc, op, x, alpha_osc, alpha_exp, n_center, n_floor, solver_type)

    class(ad_shooter_t), intent(out)         :: this
    class(mech_coeffs_t), intent(in), target :: mc
    type(oscpar_t), intent(in), target       :: op
    real(WP), intent(in)                     :: x(:)
    real(WP), intent(in)                     :: alpha_osc
    real(WP), intent(in)                     :: alpha_exp
    integer, intent(in)                      :: n_center
    integer, intent(in)                      :: n_floor
    character(LEN=*), intent(in)             :: solver_type

    ! Initialize the ad_shooter

    this%mc => mc
    this%op => op

    call this%jc%init(mc, op)
    
    this%x = x

    this%alpha_osc = alpha_osc
    this%alpha_exp = alpha_exp
    this%n_center = n_center
    this%n_floor = n_floor

    this%n = SIZE(x)
    this%n_e = this%jc%n_e

    this%solver_type = solver_type

    ! Finish

    return

  end subroutine init

!****

  subroutine shoot (this, omega, sm)

    class(ad_shooter_t), intent(in) :: this
    complex(WP), intent(in)         :: omega
    class(sysmtx_t), intent(inout)  :: sm

    integer             :: k
    complex(WP)         :: E_l(this%n_e,this%n_e)
    complex(WP)         :: E_r(this%n_e,this%n_e)
    type(ext_complex_t) :: scale

    ! Set the sysmtx equation blocks by solving IVPs across the
    ! intervals x(k) -> x(k+1)

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, scale)
    block_loop : do k = 1,this%n-1
       call solve(this%solver_type, this%jc, omega, this%x(k), this%x(k+1), E_l, E_r, scale)
       call sm%set_block(k, E_l, E_r, scale)
    end do block_loop

    ! Finish

  end subroutine shoot

!****

  subroutine recon_sh (this, omega, y_sh, x, y)

    class(ad_shooter_t), intent(in)       :: this
    complex(WP), intent(in)               :: omega
    complex(WP), intent(in)               :: y_sh(:,:)
    real(WP), intent(out), allocatable    :: x(:)
    complex(WP), intent(out), allocatable :: y(:,:)

    integer :: dn(this%n-1)
    integer :: i_a(this%n-1)
    integer :: i_b(this%n-1)
    integer :: k

    $CHECK_BOUNDS(SIZE(y_sh, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y_sh, 2),this%n)

    ! Reconstruct the eigenfunctions on a dynamically-allocated grid

    ! Allocate the grid

    dn = 0

    call plan_dispersion_grid(this%x, this%mc, omega, this%op, &
                              this%alpha_osc, this%alpha_exp, this%n_center, this%n_floor, dn)

    call build_oversamp_grid(this%x, dn, x)

    allocate(y(this%n_e,SIZE(x)))

    ! Reconstruct the eigenfunctions

    i_a(1) = 1
    i_b(1) = dn(1) + 2

    index_loop : do k = 2,this%n-1
       i_a(k) = i_b(k-1) + 1
       i_b(k) = i_a(k) + dn(k)
    end do index_loop

    !$OMP PARALLEL DO
    recon_loop : do k = 1,this%n-1
       call recon(this%solver_type, this%jc, omega, this%x(k), this%x(k+1), y_sh(:,k), y_sh(:,k+1), &
                  x(i_a(k):i_b(k)), y(:,i_a(k):i_b(k)))
    end do recon_loop
    
    ! Finish

    return

  end subroutine recon_sh

end module gyre_ad_shooter
