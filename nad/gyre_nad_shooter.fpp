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

  use gyre_base_coeffs
  use gyre_therm_coeffs
  use gyre_oscpar
  use gyre_numpar
  use gyre_ad_jacobian
  use gyre_nad_jacobian
  use gyre_sysmtx
  use gyre_ext_arith
  use gyre_ivp, ivp_abscissa => abscissa
  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: nad_shooter_t
     private
     class(base_coeffs_t), pointer  :: bc => null()
     class(therm_coeffs_t), pointer :: tc => null()
     type(oscpar_t), pointer        :: op => null()
     type(numpar_t), pointer        :: np => null()
     type(ad_jacobian_t)            :: ad_jc
     type(nad_jacobian_t)           :: nad_jc
     integer, public                :: n_e
   contains
     private
     procedure, public :: init
     procedure, public :: shoot
     procedure, public :: recon => recon_sh
     procedure, public :: abscissa
  end type nad_shooter_t

  ! Access specifiers

  private

  public :: nad_shooter_t

  ! Procedures

contains

  subroutine init (this, bc, tc, op, np)

    class(nad_shooter_t), intent(out)         :: this
    class(base_coeffs_t), intent(in), target  :: bc
    class(therm_coeffs_t), intent(in), target :: tc
    type(oscpar_t), intent(in), target        :: op
    type(numpar_t), intent(in), target        :: np

    ! Initialize the nad_shooter

    this%bc => bc
    this%tc => tc

    this%op => op
    this%np => np

    call this%ad_jc%init(bc, op)
    call this%nad_jc%init(bc, tc, op)
    
    this%n_e = this%nad_jc%n_e

    ! Finish

    return

  end subroutine init

!****

  subroutine shoot (this, omega, x, sm, x_ad)

    class(nad_shooter_t), intent(in) :: this
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x(:)
    class(sysmtx_t), intent(inout)   :: sm
    real(WP), intent(in), optional   :: x_ad

    real(WP)            :: x_ad_
    integer             :: k
    complex(WP)         :: E_l(this%n_e,this%n_e)
    complex(WP)         :: E_r(this%n_e,this%n_e)
    type(ext_complex_t) :: scale
    complex(WP)         :: lambda

    if(PRESENT(x_ad)) then
       x_ad_ = x_ad
    else
       x_ad_ = 0._WP
    endif

    ! Set the sysmtx equation blocks by solving IVPs across the
    ! intervals x(k) -> x(k+1)

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, scale, lambda) SCHEDULE (DYNAMIC)
    block_loop : do k = 1,SIZE(x)-1

       if(x(k) < x_ad_) then

          ! Shoot adiabatically

          call solve(this%np%ivp_solver_type, this%ad_jc, omega, x(k), x(k+1), E_l(1:4,1:4), E_r(1:4,1:4), scale)

          ! Fix up the thermal parts of the block

          E_l(1:4,5:6) = 0._WP
          E_r(1:4,5:6) = 0._WP

          E_l(5,:) = diff_coeffs(this%bc, this%tc, this%op, omega, x(k))
          E_r(5,:) = 0._WP
          
          E_l(6,:) = -[0._WP,0._WP,0._WP,0._WP,1._WP,0._WP]
          E_r(6,:) =  [0._WP,0._WP,0._WP,0._WP,1._WP,0._WP]

       else

          ! Shoot nonadiabatically

          call solve(this%np%ivp_solver_type, this%nad_jc, omega, x(k), x(k+1), E_l, E_r, scale)

          ! Apply the thermal-term rescaling, to assist the rootfinder

          associate(x_mid => 0.5_WP*(x(k) + x(k+1)))
            associate(V => this%bc%V(x_mid), nabla => this%tc%nabla(x_mid), &
                      c_rad => this%tc%c_rad(x_mid), c_thm => this%tc%c_thm(x_mid))
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

  subroutine recon_sh (this, omega, x_sh, y_sh, x, y, x_ad)

    class(nad_shooter_t), intent(in) :: this
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_sh(:)
    complex(WP), intent(in)          :: y_sh(:,:)
    real(WP), intent(in)             :: x(:)
    complex(WP), intent(out)         :: y(:,:)
    real(WP), intent(in), optional   :: x_ad

    real(WP)    :: x_ad_
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

    if(PRESENT(x_ad)) then
       x_ad_ = x_ad
    else
       x_ad_ = 0._WP
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

          if(x_sh(k) < x_ad_) then

             ! Reconstruct adiabatically

             call recon(this%np%ivp_solver_type, this%ad_jc, omega, x_sh(k), x_sh(k+1), y_sh(1:4,k), y_sh(1:4,k+1), &
                        x_in(:n_in), y_in(1:4,:n_in))

             y_in(5,:n_in) = y_sh(5,k)

             do i = 1,n_in
                A_5 = diff_coeffs(this%bc, this%tc, this%op, omega, x(k))
                y_in(6,i) = -DOT_PRODUCT(y_in(1:5,i), A_5(1:5))/A_5(6)
             end do
          
          else

             ! Reconstruct nonadiabatically

             call recon(this%np%ivp_solver_type, this%nad_jc, omega, x_sh(k), x_sh(k+1), y_sh(:,k), y_sh(:,k+1), &
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

    class(nad_shooter_t), intent(in) :: this
    real(WP), intent(in)             :: x_sh(:)
    real(WP), allocatable            :: x(:)

    integer :: k
    integer :: n_cell(SIZE(x_sh)-1)
    integer :: i

    ! Determine the abscissa used for shooting on the grid x_sh

    !$OMP PARALLEL DO SCHEDULE (DYNAMIC)
    count_loop : do k = 1,SIZE(x_sh)-1
       n_cell(k) = SIZE(ivp_abscissa(this%np%ivp_solver_type, x_sh(k), x_sh(k+1)))
    end do count_loop

    allocate(x(SUM(n_cell)))

    i = 1

    cell_loop : do k = 1,SIZE(x_sh)-1
       x(i:i+n_cell(k)-1) = ivp_abscissa(this%np%ivp_solver_type, x_sh(k), x_sh(k+1))
       i = i + n_cell(k)
    end do cell_loop

    $CHECK_BOUNDS(i,SIZE(x)+1)

    ! Finish

    return

  end function abscissa

!****

  function diff_coeffs (bc, tc, op, omega, x) result (a)

    class(base_coeffs_t), intent(in)  :: bc
    class(therm_coeffs_t), intent(in) :: tc
    type(oscpar_t), intent(in)        :: op
    complex(WP), intent(in)           :: omega
    real(WP), intent(in)              :: x
    complex(WP)                       :: a(6)

    ! Calculate the coefficients of the (algebraic) adiabatic
    ! diffusion equation

    associate(U => bc%U(x), c_1 => bc%c_1(x), &
              nabla_ad => bc%nabla_ad(x), &
              c_rad => tc%c_rad(x), c_dif => tc%c_dif(x), nabla => tc%nabla(x), &
              l => op%l)

      a(1) = (nabla_ad*(U - c_1*omega**2) - 4._WP*(nabla_ad - nabla) + c_dif)
      a(2) = (l*(l+1)/(c_1*omega**2)*(nabla_ad - nabla) - c_dif)
      a(3) = c_dif
      a(4) = nabla_ad
      a(5) = 0._WP
      a(6) = -nabla/c_rad
      
    end associate

    ! Finish

    return

  end function diff_coeffs

end module gyre_nad_shooter
