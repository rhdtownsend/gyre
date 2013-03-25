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

  use gyre_mech_coeffs
  use gyre_therm_coeffs
  use gyre_oscpar
  use gyre_ad_jacobian
  use gyre_nad_jacobian
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
     class(mech_coeffs_t), pointer  :: mc => null()
     class(therm_coeffs_t), pointer :: tc => null()
     type(oscpar_t), pointer        :: op => null()
     type(ad_jacobian_t)            :: ad_jc
     type(nad_jacobian_t)           :: nad_jc
     real(WP), allocatable          :: x(:)
     real(WP)                       :: alpha_osc
     real(WP)                       :: alpha_exp
     integer                        :: n_center
     integer                        :: n_floor
     integer, public                :: n
     integer, public                :: n_e
     character(LEN=256)             :: solver_type
   contains
     private
     procedure, public :: init
     procedure, public :: shoot
     procedure, public :: recon => recon_sh
  end type nad_shooter_t

  ! Access specifiers

  private

  public :: nad_shooter_t

  ! Procedures

contains

  subroutine init (this, mc, tc, op, x, alpha_osc, alpha_exp, n_center, n_floor, solver_type)

    class(nad_shooter_t), intent(out)         :: this
    class(mech_coeffs_t), intent(in), target  :: mc
    class(therm_coeffs_t), intent(in), target :: tc
    type(oscpar_t), intent(in), target        :: op
    real(WP), intent(in)                      :: x(:)
    real(WP), intent(in)                      :: alpha_osc
    real(WP), intent(in)                      :: alpha_exp
    integer, intent(in)                       :: n_center
    integer, intent(in)                       :: n_floor
    character(LEN=*), intent(in)              :: solver_type

    ! Initialize the nad_shooter

    this%mc => mc
    this%tc => tc
    this%op => op

    call this%ad_jc%init(mc, op)
    call this%nad_jc%init(mc, tc, op)
    
    this%x = x

    this%alpha_osc = alpha_osc
    this%alpha_exp = alpha_exp
    this%n_center = n_center
    this%n_floor = n_floor

    this%n = SIZE(x)
    this%n_e = this%nad_jc%n_e

    this%solver_type = solver_type

    ! Finish

    return

  end subroutine init

!****

  subroutine shoot (this, omega, sm, x_ad)

    class(nad_shooter_t), intent(in) :: this
    complex(WP), intent(in)          :: omega
    class(sysmtx_t), intent(inout)   :: sm
    real(WP), intent(in), optional   :: x_ad

    real(WP)            :: x_ad_
    integer             :: k
    complex(WP)         :: E_l(this%n_e,this%n_e)
    complex(WP)         :: E_r(this%n_e,this%n_e)
    type(ext_complex_t) :: scale
    complex(WP)         :: A(this%n_e,this%n_e)
    complex(WP)         :: lambda

    if(PRESENT(x_ad)) then
       x_ad_ = x_ad
    else
       x_ad_ = 0._WP
    endif

    ! Set the sysmtx equation blocks by solving IVPs across the
    ! intervals x(k) -> x(k+1)

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, scale, A, lambda)
    block_loop : do k = 1,this%n-1

       if(this%x(k) < x_ad_) then

          ! Shoot adiabatically

          call solve(this%solver_type, this%ad_jc, omega, this%x(k), this%x(k+1), E_l(1:4,1:4), E_r(1:4,1:4), scale)

          ! Fix up the thermal parts of the block

          call this%nad_jc%eval_logx(omega, this%x(k), A)

          E_l(1:4,5:6) = 0._WP
          E_r(1:4,5:6) = 0._WP

          E_l(5,:) = A(5,:)
          E_r(5,:) = 0._WP
          
          E_l(6,:) = -[0._WP,0._WP,0._WP,0._WP,1._WP,0._WP]
          E_r(6,:) =  [0._WP,0._WP,0._WP,0._WP,1._WP,0._WP]

       else

          ! Shoot nonadiabatically

          call solve(this%solver_type, this%nad_jc, omega, this%x(k), this%x(k+1), E_l, E_r, scale)

          ! Apply the thermal-term rescaling, to assist the rootfinder

          associate(x_mid => 0.5_WP*(this%x(k) + this%x(k+1)))
            associate(V => this%mc%V(x_mid), nabla => this%tc%nabla(x_mid), &
                      c_rad => this%tc%c_rad(x_mid), c_thm => this%tc%c_thm(x_mid))
              lambda = SQRT(V*nabla/c_rad * (0._WP,1._WP)*omega*c_thm)/x_mid
            end associate
          end associate

          scale = scale*exp(ext_complex(-lambda*(this%x(k+1)-this%x(k))))

       endif

       call sm%set_block(k, E_l, E_r, scale)

    end do block_loop

    ! Finish

  end subroutine shoot

!****

  subroutine recon_sh (this, omega, y_sh, x, y, x_ad)

    class(nad_shooter_t), intent(in)      :: this
    complex(WP), intent(in)               :: omega
    complex(WP), intent(in)               :: y_sh(:,:)
    real(WP), intent(out), allocatable    :: x(:)
    complex(WP), intent(out), allocatable :: y(:,:)
    real(WP), intent(in), optional        :: x_ad

    real(WP)    :: x_ad_
    integer     :: dn(this%n-1)
    integer     :: i_a(this%n-1)
    integer     :: i_b(this%n-1)
    integer     :: k
    integer     :: i
    complex(WP) :: A(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(y_sh, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y_sh, 2),this%n)

    if(PRESENT(x_ad)) then
       x_ad_ = x_ad
    else
       x_ad_ = 0._WP
    endif

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

       if(this%x(k) < x_ad_) then

          ! Reconstruct adiabatically

          call recon(this%solver_type, this%ad_jc, omega, this%x(k), this%x(k+1), y_sh(1:4,k), y_sh(1:4,k+1), &
                     x(i_a(k):i_b(k)), y(1:4,i_a(k):i_b(k)))

          y(5,i_a(k):i_b(k)) = y_sh(5,k)

          do i = i_a(k),i_b(k)
             call this%nad_jc%eval_logx(omega, this%x(i), A)
             y(6,i) = -DOT_PRODUCT(y(1:5,i), A(5,1:5))/A(5,6)
          end do
          
       else

          ! Reconstruct nonadiabatically

          call recon(this%solver_type, this%nad_jc, omega, this%x(k), this%x(k+1), y_sh(:,k), y_sh(:,k+1), &
                     x(i_a(k):i_b(k)), y(:,i_a(k):i_b(k)))

       endif

    end do recon_loop
    
    ! Finish

    return

  end subroutine recon_sh

end module gyre_nad_shooter
