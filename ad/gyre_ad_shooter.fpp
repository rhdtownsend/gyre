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
  use gyre_numpar
  use gyre_ad_jacobian
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_ivp

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ad_shooter_t
     private
     class(mech_coeffs_t), pointer :: mc => null()
     type(oscpar_t), pointer       :: op => null()
     type(numpar_t), pointer       :: np => null()
     type(ad_jacobian_t)           :: jc
     integer, public               :: n_e
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

  subroutine init (this, mc, op, np)

    class(ad_shooter_t), intent(out)         :: this
    class(mech_coeffs_t), intent(in), target :: mc
    type(oscpar_t), intent(in), target       :: op
    type(numpar_t), intent(in), target       :: np

    ! Initialize the ad_shooter

    this%mc => mc
    this%op => op
    this%np => np

    call this%jc%init(mc, op)
    
    this%n_e = this%jc%n_e

    ! Finish

    return

  end subroutine init

!****

  subroutine shoot (this, omega, x, sm)

    class(ad_shooter_t), intent(in) :: this
    complex(WP), intent(in)         :: omega
    real(WP), intent(in)            :: x(:)
    class(sysmtx_t), intent(inout)  :: sm

    integer             :: k
    complex(WP)         :: E_l(this%n_e,this%n_e)
    complex(WP)         :: E_r(this%n_e,this%n_e)
    type(ext_complex_t) :: S

    ! Set the sysmtx equation blocks by solving IVPs across the
    ! intervals x(k) -> x(k+1)

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, S)
    block_loop : do k = 1,SIZE(x)-1
       call solve(this%np%ivp_solver_type, this%jc, omega, x(k), x(k+1), E_l, E_r, S)
       call sm%set_block(k, E_l, E_r, S)
    end do block_loop

    ! Finish

  end subroutine shoot

!****

  subroutine recon_sh (this, omega, x_sh, y_sh, x, y)

    class(ad_shooter_t), intent(in) :: this
    complex(WP), intent(in)         :: omega
    real(WP), intent(in)            :: x_sh(:)
    complex(WP), intent(in)         :: y_sh(:,:)
    real(WP), intent(in)            :: x(:)
    complex(WP), intent(out)        :: y(:,:)

    integer     :: n_sh
    integer     :: n
    integer     :: k
    logical     :: mask(SIZE(x))
    integer     :: n_in
    integer     :: i
    integer     :: i_in(SIZE(x))
    real(WP)    :: x_in(SIZE(x))
    complex(WP) :: y_in(this%n_e,SIZE(x))

    $CHECK_BOUNDS(SIZE(y_sh, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y_sh, 2),SIZE(x_sh))

    $CHECK_BOUNDS(SIZE(y, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Reconstruct the eigenfunctions on the supplied grid

    n_sh = SIZE(x_sh)
    n = SIZE(x)

    !$OMP PARALLEL DO PRIVATE (mask, n_in, i_in, x_in, y_in)
    recon_loop : do k = 1,n_sh-1

       ! Select those points which fall in the current interval

       if(k == 1) then
          mask = x < x_sh(k+1)
       elseif(k == n_sh-1) then
          mask = x >= x_sh(k)
       else
          mask = x >= x_sh(k) .AND. x < x_sh(k+1)
       endif

       n_in = COUNT(mask)

       if(n_in > 0) then

          ! Reconstruct in the interval

          i_in(:n_in) = PACK([(i,i=1,n)], MASK=mask)

          x_in(:n_in) = x(i_in(:n_in))

          call recon(this%np%ivp_solver_type, this%jc, omega, x_sh(k), x_sh(k+1), y_sh(:,k), y_sh(:,k+1), &
               x_in(:n_in), y_in(:,:n_in))

          y(:,i_in(:n_in)) = y_in(:,:n_in)

       endif

    end do recon_loop
    
    ! Finish

    return

  end subroutine recon_sh

end module gyre_ad_shooter
