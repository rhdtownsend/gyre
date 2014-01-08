! Module   : gyre_shooter_ad
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

module gyre_shooter_ad

  ! Uses

  use core_kinds

  use gyre_coeffs
  use gyre_oscpar
  use gyre_numpar
  use gyre_ivp
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_ivp

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: shooter_ad_t
     private
     class(coeffs_t), pointer  :: cf => null()
     class(ivp_t), allocatable :: iv
     type(oscpar_t), pointer   :: op
     type(numpar_t), pointer   :: np
     integer, public           :: n_e
   contains
     private
     procedure, public :: shoot
     procedure, public :: recon => recon_sh
     procedure, public :: abscissa
  end type shooter_ad_t

  ! Interfaces

  interface shooter_ad_t
     module procedure init_sh
  end interface shooter_ad_t

  ! Access specifiers

  private

  public :: shooter_ad_t

  ! Procedures

contains

  function init_sh (cf, iv, op, np) result (sh)

    class(coeffs_t), pointer, intent(in) :: cf
    class(ivp_t), intent(in)             :: iv
    type(oscpar_t), intent(in)           :: op
    type(numpar_t), intent(in)           :: np
    type(shooter_ad_t)                   :: sh

    ! Construct the shooter_ad

    sh%cf => cf
    allocate(sh%iv, SOURCE=iv)
    sh%op = op
    sh%np = np

    sh%n_e = sh%iv%n_e

    ! Finish

    return

  end function init_sh

!****

  subroutine shoot (this, omega, x, sm)

    class(shooter_ad_t), intent(in) :: this
    complex(WP), intent(in)         :: omega
    real(WP), intent(in)            :: x(:)
    class(sysmtx_t), intent(inout)  :: sm

    integer             :: k
    complex(WP)         :: E_l(this%n_e,this%n_e)
    complex(WP)         :: E_r(this%n_e,this%n_e)
    type(ext_complex_t) :: S

    ! Set the sysmtx equation blocks by solving IVPs across the
    ! intervals x(k) -> x(k+1)

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, S) SCHEDULE (DYNAMIC)
    block_loop : do k = 1,SIZE(x)-1
       call this%iv%solve(omega, x(k), x(k+1), E_l, E_r, S, use_real=.TRUE.)
       call sm%set_block(k, E_l, E_r, S)
    end do block_loop

    ! Finish

  end subroutine shoot

!****

  subroutine recon_sh (this, omega, x_sh, y_sh, x, y)

    class(shooter_ad_t), intent(in) :: this
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

    !$OMP PARALLEL DO PRIVATE (mask, n_in, i_in, x_in, y_in) SCHEDULE (DYNAMIC)
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

          call this%iv%recon(omega, x_sh(k), x_sh(k+1), y_sh(:,k), y_sh(:,k+1), &
               x_in(:n_in), y_in(:,:n_in), use_real=.TRUE.)

          y(:,i_in(:n_in)) = y_in(:,:n_in)

       endif

    end do recon_loop
    
    ! Finish

    return

  end subroutine recon_sh

!****

  function abscissa (this, x_sh) result (x)

    class(shooter_ad_t), intent(in) :: this
    real(WP), intent(in)            :: x_sh(:)
    real(WP), allocatable           :: x(:)

    integer :: k
    integer :: n_cell(SIZE(x_sh)-1)
    integer :: i

    ! Determine the abscissa used for shooting on the grid x_sh

    !$OMP PARALLEL DO SCHEDULE (DYNAMIC)
    count_loop : do k = 1,SIZE(x_sh)-1
       n_cell(k) = SIZE(this%iv%abscissa(x_sh(k), x_sh(k+1)))
    end do count_loop

    allocate(x(SUM(n_cell)))

    i = 1

    cell_loop : do k = 1,SIZE(x_sh)-1
       x(i:i+n_cell(k)-1) = this%iv%abscissa(x_sh(k), x_sh(k+1))
       i = i + n_cell(k)
    end do cell_loop

    $CHECK_BOUNDS(i,SIZE(x)+1)

    ! Finish

    return

  end function abscissa

end module gyre_shooter_ad
