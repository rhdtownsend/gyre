! Module   : gyre_grid
! Purpose  : grid construction
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

module gyre_grid

  ! Uses

  use core_kinds
  use core_constants

  use gyre_mech_coeffs
  use gyre_oscpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: build_geom_grid
  public :: build_log_grid
  public :: build_oversamp_grid
  public :: plan_dispersion_grid

  ! Procedures

contains

  subroutine build_geom_grid (s, n, x)

    real(WP), intent(in)                :: s
    integer, intent(in)                 :: n
    real(WP), intent(out), allocatable  :: x(:)

    integer  :: m
    integer  :: k
    real(WP) :: dx_1

    ! Build an n-point grid with geometric spacing in each half of the
    ! [0,1] interval. The parameter s controls the ratio between the
    ! central cell size and the boundary cell size

    allocate(x(n))

    if(MOD(n, 2) == 0) then

       ! Even number of grid points / odd number of cells

       ! Solve for the size of the boundary cells

       m = n/2

       $if($GFORTRAN_PR_56872)
       dx_1 = 0.5_WP/(SUM(s**[(REAL(k-1, WP)/REAL(m-1, WP),k=1,m)]) - &
                      0.5_WP*s)
       $else
       dx_1 = 0.5_WP/(SUM([(s**(REAL(k-1, WP)/REAL(m-1, WP)),k=1,m)]) - &
                      0.5_WP*s)
       $endif

       ! Set up the inner part of the grid

       x(1) = 0._WP
       
       even_grid_loop : do k = 2,m
          x(k) = x(k-1) + dx_1*s**(REAL(k-2, WP)/REAL(m-1, WP))
       end do even_grid_loop

       ! Reflect to get the outer part of the grid

       x(m+1:) = 1._WP - x(m:1:-1)

    else

       ! Odd number of grid points / even number of cells

       ! Solve for the size of the boundary cells

       m = (n-1)/2

       $if($GFORTRAN_PR_56872)
       dx_1 = 0.5_WP/(SUM(s**[(REAL(k-1, WP)/REAL(m-1, WP),k=1,m)]))
       $else
       dx_1 = 0.5_WP/(SUM([(s**(REAL(k-1, WP)/REAL(m-1, WP)),k=1,m)]))
       $endif

       ! Set up the inner part of the grid

       x(1) = 0._WP
       
       odd_grid_loop : do k = 2,m
          x(k) = x(k-1) + dx_1*s**(REAL(k-2, WP)/REAL(m-1, WP))
       end do odd_grid_loop

       x(m+1) = 0.5_WP

       ! Reflect to get the outer part of the grid

       x(m+2:) = 1._WP - x(m:1:-1)

    end if

    ! Finish

    return

  end subroutine build_geom_grid

!****

  subroutine build_log_grid (s, n, x)

    real(WP), intent(in)                :: s
    integer, intent(in)                 :: n
    real(WP), intent(out), allocatable  :: x(:)

    real(WP) :: dx_1
    integer  :: k
    real(WP) :: w
    real(WP) :: t

    ! Build an n-point grid with logarithmic spacing in each half of
    ! the [0,1] interval. The parameter s controls the ratio between
    ! the mean cell size and the boundary cell size. omega is not used

    allocate(x(n))

    dx_1 = 1._WP/(s*(n-1))

    if(MOD(n, 2) == 0) then

       ! Even number of grid points / odd number of cells

       ! Set up the inner part of the grid

       x(1) = 0._WP

       even_grid_loop : do k = 2,n/2

          w = (k-1.5_WP)/(n/2-1.5_WP)
          t = LOG(0.5_WP)*(1._WP-w) + LOG(dx_1)*w

          x(n/2-k+2) = EXP(t)

       enddo even_grid_loop
       
       ! Reflect to get the outer part of the grid

       x(n/2+1:) = 1._WP - x(n/2:1:-1)

    else

       ! Odd number of grid points / even number of cells

       ! Set up the inner part of the grid

       x(1) = 0._WP

       odd_grid_loop : do k = 2,(n-1)/2

          w = (k-1._WP)/((n-1)/2-1._WP)
          t = LOG(0.5_WP)*(1._WP-w) + LOG(dx_1)*w

          x((n-1)/2-k+2) = EXP(t)

       end do odd_grid_loop

       x((n+1)/2) = 0.5_WP

       ! Reflect to get the outer part of the grid

       x((n+1)/2+1:) = 1._WP - x((n-1)/2:1:-1)

    end if

    ! Finish

    return

  end subroutine build_log_grid

!****

  subroutine build_oversamp_grid (x_in, dn, x)

    real(WP), intent(in)               :: x_in(:)
    integer, intent(in)                :: dn(:)
    real(WP), allocatable, intent(out) :: x(:)

    integer :: n
    integer :: i
    integer :: j
    integer :: k
    
    $CHECK_BOUNDS(SIZE(dn),SIZE(x_in)-1)

    ! Build a grid by oversampling x_in, with dn(i) additional points
    ! inserted uniformly across cell i

    n = SIZE(x_in)

    allocate(x(SUM(dn) + n))

    k = 1

    do i = 1,n-1
       do j = 1,dn(i)+1
          x(k) = x_in(i) + (j-1)*(x_in(i+1)-x_in(i))/(dn(i)+1)
          k = k + 1
       end do
    end do
    
    x(k) = x_in(n)

    ! Finish

    return

  end subroutine build_oversamp_grid

!****

  subroutine plan_dispersion_grid (x_mc, mc, omega, op, alpha_osc, alpha_exp, n_center, n_floor, dn)

    real(WP), intent(in)             :: x_mc(:)
    class(mech_coeffs_t), intent(in) :: mc
    complex(WP), intent(in)          :: omega
    type(oscpar_t), intent(in)       :: op
    real(WP), intent(in)             :: alpha_osc
    real(WP), intent(in)             :: alpha_exp
    integer, intent(in)              :: n_center
    integer, intent(in)              :: n_floor
    integer, intent(inout)           :: dn(:)

    integer     :: n
    integer     :: i
    complex(WP) :: b
    complex(WP) :: c
    complex(WP) :: k_r_pos(SIZE(x_mc))
    complex(WP) :: k_r_neg(SIZE(x_mc))
    real(WP)    :: dphi_osc
    real(WP)    :: dphi_exp
    integer     :: i_turn
    real(WP)    :: x_turn

    $CHECK_BOUNDS(SIZE(dn),SIZE(x_mc)-1)

    ! Plan an oversamp grid, modifying dn (see build_grid_oversamp)
    ! with additional points based on a local dispersion analysis

    n = SIZE(x_mc)

    ! Estimate local radial wavenumbers (note that only the real part
    ! of omega is used)

    k_r_pos(1) = 0._WP
    k_r_neg(1) = 0._WP

    $if(!$GFORTRAN_PR_56052)
    !$OMP PARALLEL DO PRIVATE (b, c)
    $endif
    wave_loop : do i = 2,n

       associate(x => x_mc(i))
         associate(V_g => mc%V(x)/mc%Gamma_1(x), As => mc%As(x), U => mc%U(x), c_1 => mc%c_1(x), &
                   l => op%l, omega_re => REAL(omega))

           b = CMPLX(0._WP, 1._WP, WP)*(V_g + As - U - 2._WP - 2._WP*(l-2))
           c = (l*(l+1)/(c_1*omega_re**2) - V_g)*(c_1*omega_re**2 - As)

           k_r_pos(i) = (0.5_WP*(-b + SQRT(b**2 - 4._WP*c)))/x
           k_r_neg(i) = (0.5_WP*(-b - SQRT(b**2 - 4._WP*c)))/x

         end associate
       end associate

    end do wave_loop

    ! Place points to ensure a given sampling of
    ! oscillatory/exponential scale lengths

    !$OMP PARALLEL DO PRIVATE (dphi_osc, dphi_exp)
    samp_loop : do i = 1,n-1

       dphi_osc = MAX(ABS(REAL(k_r_pos(i))), ABS(REAL(k_r_pos(i+1))), &
                      ABS(REAL(k_r_neg(i))), ABS(REAL(k_r_neg(i+1))))*(x_mc(i+1) - x_mc(i))
       dphi_exp = MAX(ABS(AIMAG(k_r_pos(i))), ABS(AIMAG(k_r_pos(i+1))), &
                      ABS(AIMAG(k_r_neg(i))), ABS(AIMAG(k_r_neg(i+1))))*(x_mc(i+1) - x_mc(i))

       dn(i) = MAX(dn(i), FLOOR((alpha_osc*dphi_osc)/PI), FLOOR((alpha_exp*dphi_exp)/PI))

    end do samp_loop

    ! Place points to ensure the central evanescent zone has at least
    ! n_center points in it

    ! First, locate the innermost turning point

    i_turn = 1

    turn_loop : do i = 2,n
       if(REAL(k_r_pos(i)) == 0._WP .AND. REAL(k_r_neg(i)) == 0._WP) then
          i_turn = i
       else
          exit turn_loop
       endif
    end do turn_loop

    ! Place points

    if(i_turn == 1) then

       ! Turning point is in innermost cell; assume V_g, As ~ x^2 to
       ! interpolate its position

       associate(x => x_mc(2))

         associate(V_g => mc%V(x)/mc%Gamma_1(x), As => mc%As(x), c_1 => mc%c_1(x), &
                   l => op%l, omega_re => REAL(omega))
           if(l /= 0._WP) then
              if(As > 0._WP) then
                 x_turn = SQRT(MIN(c_1*omega_re**2/As, l*(l+1)/(c_1*omega_re**2*V_g)))*x
              else
                 x_turn = SQRT(l*(l+1)/(c_1*omega_re**2*V_g))*x
              endif
           else
              if(As > 0._WP) then
                 x_turn = SQRT(c_1*omega_re**2/As)
              else
                 x_turn = HUGE(0._WP)
              endif
           endif
         end associate

         dn(1) = MAX(dn(1), CEILING(x/x_turn*n_center))

       end associate

    else

       dn(:i_turn-1) = MAX(dn(:i_turn-1), CEILING(REAL(n_center, WP)/(i_turn-1)))

    endif

    ! Place points based on a simple floor

    dn = MAX(dn, n_floor)

    ! Finish

    return

  end subroutine plan_dispersion_grid

end module gyre_grid
