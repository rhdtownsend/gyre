! Module   : gyre_shooter_nad
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

module gyre_shooter_nad

  ! Uses

  use core_kinds
  use core_order

  use gyre_coeffs
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

  type :: shooter_nad_t
     private
     class(coeffs_t), pointer  :: cf => null()
     class(ivp_t), allocatable :: iv
     class(ivp_t), allocatable :: iv_upw
     type(oscpar_t)            :: op
     type(numpar_t)            :: np
     integer, public           :: n_e
   contains
     private
     procedure, public :: shoot
     procedure, public :: recon => recon_sh
     procedure, public :: abscissa
  end type shooter_nad_t

  ! Interfaces

  interface shooter_nad_t
     module procedure init_sh
  end interface shooter_nad_t

  ! Access specifiers

  private

  public :: shooter_nad_t

  ! Procedures

contains

  function init_sh (cf, iv, iv_upw, op, np) result (sh)

    class(coeffs_t), pointer, intent(in) :: cf
    class(ivp_t), intent(in)             :: iv
    class(ivp_t), intent(in)             :: iv_upw
    type(oscpar_t), intent(in)           :: op
    type(numpar_t), intent(in)           :: np
    type(shooter_nad_t)                  :: sh

    ! Construct the shooter_nad

    sh%cf => cf
    allocate(sh%iv, SOURCE=iv)
    allocate(sh%iv_upw, SOURCE=iv_upw)
    sh%op = op
    sh%np = np

    sh%n_e = sh%iv%n_e

    ! Finish

    return

  end function init_sh

!****

  subroutine shoot (this, omega, x, sm, x_upw)

    class(shooter_nad_t), intent(in) :: this
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x(:)
    class(sysmtx_t), intent(inout)   :: sm
    real(WP), intent(in), optional   :: x_upw

    real(WP)            :: x_upw_
    integer             :: k
    complex(WP)         :: E_l(this%n_e,this%n_e)
    complex(WP)         :: E_r(this%n_e,this%n_e)
    type(ext_complex_t) :: scale
    complex(WP)         :: lambda

    if(PRESENT(x_upw)) then
       x_upw_ = x_upw
    else
       x_upw_ = 0._WP
    endif

    ! Set the sysmtx equation blocks by solving IVPs across the
    ! intervals x(k) -> x(k+1)

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, scale, lambda) SCHEDULE (DYNAMIC)
    block_loop : do k = 1,SIZE(x)-1

       if(k > 1 .AND. x(k) < x_upw_) then

          ! Shoot using upwinded finite differences

          call this%iv_upw%solve(omega, x(k), x(k+1), E_l, E_r, scale)

       else

          ! Shoot

          call this%iv%solve(omega, x(k), x(k+1), E_l, E_r, scale)

          ! Apply the thermal-term rescaling, to assist the rootfinder

          associate(x_mid => 0.5_WP*(x(k) + x(k+1)))
            associate(V => this%cf%V(x_mid), nabla => this%cf%nabla(x_mid), &
                      c_rad => this%cf%c_rad(x_mid), c_thm => this%cf%c_thm(x_mid))
              lambda = SQRT(V*nabla/c_rad * (0._WP,1._WP)*omega*c_thm)/x_mid
            end associate
          end associate

          scale = scale*exp(ext_complex(-lambda*(x(k+1)-x(k))))

       endif

       call sm%set_block(k, E_l, E_r, scale)

    end do block_loop

    ! Finish

  end subroutine shoot

!****

  subroutine recon_sh (this, omega, x_sh, y_sh, x, y, x_upw)

    class(shooter_nad_t), intent(in) :: this
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_sh(:)
    complex(WP), intent(in)          :: y_sh(:,:)
    real(WP), intent(in)             :: x(:)
    complex(WP), intent(out)         :: y(:,:)
    real(WP), intent(in), optional   :: x_upw

    real(WP)    :: x_upw_
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

    if(PRESENT(x_upw)) then
       x_upw_ = x_upw
    else
       x_upw_ = 0._WP
    endif

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

          if(k > 1 .AND. x_sh(k) < x_upw_) then

             ! Reconstruct using upwinded finite differences

             call this%iv_upw%recon(omega, x_sh(k), x_sh(k+1), y_sh(:,k), y_sh(:,k+1), &
                                    x_in(:n_in), y_in(:,:n_in))

          else

             ! Reconstruct

             call this%iv%recon(omega, x_sh(k), x_sh(k+1), y_sh(:,k), y_sh(:,k+1), &
                                x_in(:n_in), y_in(:,:n_in))

          endif

          y(:,i_in(:n_in)) = y_in(:,:n_in)

       end if

    end do recon_loop

    ! Finish

    return

  end subroutine recon_sh

!****

  function abscissa (this, x_sh) result (x)

    class(shooter_nad_t), intent(in) :: this
    real(WP), intent(in)             :: x_sh(:)
    real(WP), allocatable            :: x(:)

    integer               :: k
    integer               :: n_cell(SIZE(x_sh)-1)
    integer               :: n_cell_upw(SIZE(x_sh)-1)
    real(WP), allocatable :: x_(:)
    integer               :: i

    ! Determine the abscissa used for shooting on the grid x_sh

    !$OMP PARALLEL DO SCHEDULE (DYNAMIC)
    count_loop : do k = 1,SIZE(x_sh)-1
       n_cell(k) = SIZE(this%iv%abscissa(x_sh(k), x_sh(k+1)))
       n_cell_upw(k) = SIZE(this%iv_upw%abscissa(x_sh(k), x_sh(k+1)))
    end do count_loop

    allocate(x_(SUM(n_cell+n_cell_upw)))

    i = 1

    cell_loop : do k = 1,SIZE(x_sh)-1

       x_(i:i+n_cell(k)-1) = this%iv%abscissa(x_sh(k), x_sh(k+1))
       x_(i+n_cell(k):i+n_cell_upw(k)-1) = this%iv_upw%abscissa(x_sh(k), x_sh(k+1))

       i = i + n_cell(k) + n_cell_upw(k)

    end do cell_loop

    $CHECK_BOUNDS(i,SIZE(x_)+1)

    x = x_(unique_indices(x_))

    ! Finish

    return

  end function abscissa

! !****

!   function diff_coeffs (cf, op, omega, x) result (a)

!     class(coeffs_t), intent(in) :: cf
!     type(oscpar_t), intent(in)  :: op
!     complex(WP), intent(in)     :: omega
!     real(WP), intent(in)        :: x
!     complex(WP)                 :: a(6)

!     ! Calculate the coefficients of the (algebraic) adiabatic
!     ! diffusion equation

!     associate(U => cf%U(x), c_1 => cf%c_1(x), &
!               nabla_ad => cf%nabla_ad(x), &
!               c_rad => cf%c_rad(x), c_dif => cf%c_dif(x), nabla => cf%nabla(x), &
!               l => op%l)

!       a(1) = (nabla_ad*(U - c_1*omega**2) - 4._WP*(nabla_ad - nabla) + c_dif)
!       a(2) = (l*(l+1)/(c_1*omega**2)*(nabla_ad - nabla) - c_dif)
!       a(3) = c_dif
!       a(4) = nabla_ad
!       a(5) = 0._WP
!       a(6) = -nabla/c_rad
      
!     end associate

!     ! Finish

!     return

!   end function diff_coeffs

end module gyre_shooter_nad
