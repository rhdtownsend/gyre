! Module   : gyre_nad_shooter
! Purpose  : nonadiabatic multiple shooting
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

module gyre_nad_shooter

  ! Uses

  use core_kinds
  use core_order

  use gyre_model
  use gyre_oscpar
  use gyre_numpar
  use gyre_jacobian
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_ivp
  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: nad_shooter_t
     private
     class(model_t), pointer   :: ml => null()
     class(ivp_t), allocatable :: iv
     type(numpar_t)            :: np
     integer, public           :: n_e
   contains
     private
     procedure, public :: shoot => shoot_
     procedure, public :: recon => recon_
     procedure, public :: abscissa => abscissa_
  end type nad_shooter_t

  ! Interfaces

  interface nad_shooter_t
     module procedure nad_shooter_t_
  end interface nad_shooter_t

  ! Access specifiers

  private

  public :: nad_shooter_t

  ! Procedures

contains

  function nad_shooter_t_ (ml, iv, np) result (sh)

    class(model_t), pointer, intent(in) :: ml
    class(ivp_t), intent(in)            :: iv
    type(numpar_t), intent(in)          :: np
    type(nad_shooter_t)                 :: sh

    ! Construct the nad_shooter_t

    sh%ml => ml
    allocate(sh%iv, SOURCE=iv)
    sh%np = np

    sh%n_e = sh%iv%n_e

    ! Finish

    return

  end function nad_shooter_t_

!****

  subroutine shoot_ (this, omega, x, sm)

    class(nad_shooter_t), intent(in) :: this
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x(:)
    class(sysmtx_t), intent(inout)   :: sm

    integer             :: k
    complex(WP)         :: E_l(this%n_e,this%n_e)
    complex(WP)         :: E_r(this%n_e,this%n_e)
    type(ext_complex_t) :: scale
    complex(WP)         :: lambda

    ! Set the sysmtx equation blocks by solving IVPs across the
    ! intervals x(k) -> x(k+1)

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, scale, lambda) SCHEDULE (DYNAMIC)
    block_loop : do k = 1,SIZE(x)-1

       ! Shoot

       call this%iv%solve(omega, x(k), x(k+1), E_l, E_r, scale)

       ! Apply the thermal-term rescaling, to assist the rootfinder

       associate(x_mid => x(k) + 0.5_WP*(x(k+1) - x(k)))
         associate(V => this%ml%V(x_mid), nabla => this%ml%nabla(x_mid), &
                   c_rad => this%ml%c_rad(x_mid), c_thm => this%ml%c_thm(x_mid))
           lambda = SQRT(V*nabla/c_rad * (0._WP,1._WP)*omega*c_thm)/x_mid
         end associate
       end associate

       scale = scale*exp(ext_complex_t(-lambda*(x(k+1)-x(k))))

       call sm%set_block(k, E_l, E_r, scale)

    end do block_loop

    ! Finish

  end subroutine shoot_

!****

  subroutine recon_ (this, omega, x_sh, y_sh, x, y)

    class(nad_shooter_t), intent(in) :: this
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_sh(:)
    complex(WP), intent(in)          :: y_sh(:,:)
    real(WP), intent(in)             :: x(:)
    complex(WP), intent(out)         :: y(:,:)

    integer     :: n_sh
    integer     :: n
    integer     :: k
    logical     :: mask(SIZE(x))
    integer     :: n_in
    integer     :: i
    integer     :: i_in(SIZE(x))
    real(WP)    :: x_in(SIZE(x))
    complex(WP) :: y_in(this%n_e,SIZE(x))
    complex(WP) :: A_5(this%n_e)

    $CHECK_BOUNDS(SIZE(y_sh, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y_sh, 2),SIZE(x_sh))

    $CHECK_BOUNDS(SIZE(y, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Reconstruct the eigenfunctions on the supplied grid

    n_sh = SIZE(x_sh)
    n = SIZE(x)

    !$OMP PARALLEL DO PRIVATE (mask, n_in, i_in, x_in, y_in, A_5) SCHEDULE (DYNAMIC)
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
               x_in(:n_in), y_in(:,:n_in))

          y(:,i_in(:n_in)) = y_in(:,:n_in)

       end if

    end do recon_loop

    ! Finish

    return

  end subroutine recon_

!****

  function abscissa_ (this, x_sh) result (x)

    class(nad_shooter_t), intent(in) :: this
    real(WP), intent(in)             :: x_sh(:)
    real(WP), allocatable            :: x(:)

    integer :: k
    integer :: n_cell(SIZE(x_sh)-1)
    integer :: i

    ! Determine the abscissa used for shooting on the grid x_sh

    !$OMP PARALLEL DO SCHEDULE (DYNAMIC)
    count_loop : do k = 1,SIZE(x_sh)-1
       n_cell(k) = SIZE(this%iv%abscissa(x_sh(k), x_sh(k+1)))
    end do count_loop

    allocate(x(SUM(n_cell+1)))

    i = 1

    cell_loop : do k = 1,SIZE(x_sh)-1
       x(i:i+n_cell(k)-1) = this%iv%abscissa(x_sh(k), x_sh(k+1))
       x(i+n_cell(k)) = x_sh(k) + 0.5_WP*(x_sh(k+1) - x_sh(k))
       i = i + n_cell(k) + 1
    end do cell_loop

    $CHECK_BOUNDS(i,SIZE(x)+1)

    x = x(unique_indices(x))

    ! Finish

    return

  end function abscissa_

end module gyre_nad_shooter
