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
  use core_func
  use core_order

  use gyre_base_coeffs
  use gyre_oscpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (func_t) :: gamma_t
     class(base_coeffs_t), pointer :: bc => null()
     type(oscpar_t), pointer       :: op => null()
     real(WP)                      :: omega
   contains
     procedure :: eval_c => eval_gamma
  end type gamma_t

  ! Access specifiers

  private

  public :: build_geom_grid
  public :: build_log_grid
  public :: build_dispersion_grid

  ! Procedures

contains

  subroutine build_geom_grid (s, n, x)

    real(WP), intent(in)               :: s
    integer, intent(in)                :: n
    real(WP), allocatable, intent(out) :: x(:)

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

    real(WP), intent(in)               :: s
    integer, intent(in)                :: n
    real(WP), allocatable, intent(out) :: x(:)

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

  subroutine build_dispersion_grid (x_bc, bc, op, omega_a, omega_b, alpha_osc, alpha_exp, n_center, n_floor, x)

    real(WP), intent(in)               :: x_bc(:)
    class(base_coeffs_t), intent(in)   :: bc
    type(oscpar_t), intent(in)         :: op
    real(WP), intent(in)               :: omega_a
    real(WP), intent(in)               :: omega_b
    real(WP), intent(in)               :: alpha_osc
    real(WP), intent(in)               :: alpha_exp
    integer, intent(in)                :: n_center
    integer, intent(in)                :: n_floor
    real(WP), allocatable, intent(out) :: x(:)

    integer :: dn(SIZE(x_bc)-1)

    ! Plan the grid

    call plan_dispersion_grid (x_bc, bc, op, omega_a, omega_b, alpha_osc, alpha_exp, n_center, n_floor, dn)

    ! Build it

    call build_oversamp_grid(x_bc, dn, x)

    ! Finish

    return

  end subroutine build_dispersion_grid

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

  subroutine plan_dispersion_grid (x_bc, bc, op, omega_a, omega_b, alpha_osc, alpha_exp, n_center, n_floor, dn)

    real(WP), intent(in)             :: x_bc(:)
    class(base_coeffs_t), intent(in) :: bc
    type(oscpar_t), intent(in)       :: op
    real(WP), intent(in)             :: omega_a
    real(WP), intent(in)             :: omega_b
    real(WP), intent(in)             :: alpha_osc
    real(WP), intent(in)             :: alpha_exp
    integer, intent(in)              :: n_center
    integer, intent(in)              :: n_floor
    integer, intent(out)             :: dn(:)

    integer       :: n
    real(WP)      :: lambda_osc(SIZE(x_bc))
    real(WP)      :: lambda_exp(SIZE(x_bc))
    integer       :: i
    real(WP)      :: g_4
    real(WP)      :: g_2
    real(WP)      :: g_0
    real(WP)      :: omega2_ext
    real(WP)      :: gamma_ext
    real(WP)      :: gamma_a
    real(WP)      :: gamma_b
    real(WP)      :: dphi_osc
    real(WP)      :: dphi_exp
    real(WP)      :: x_turn_a
    real(WP)      :: x_turn_b
    integer       :: i_turn
    real(WP)      :: x_turn

    $CHECK_BOUNDS(SIZE(dn),SIZE(x_bc)-1)

    $ASSERT(omega_b >= omega_a,Incorrect frequency ordering)

    ! Plan an oversamp grid, modifying dn (see build_grid_oversamp)
    ! with additional points based on a local dispersion analysis

    dn = 0

    n = SIZE(x_bc)

    ! Calculate the maximal oscillatory and exponential wavenumbers at
    ! each point of x_bc for (real) frequencies between omega_a and
    ! omega_b

    lambda_osc(1) = 0._WP
    lambda_exp(1) = 0._WP

    wave_loop : do i = 2,n-1

       associate(x => x_bc(i))
         associate(V_g => bc%V(x)/bc%Gamma_1(x), As => bc%As(x), U => bc%U(x), c_1 => bc%c_1(x), &
                   l => op%l)

           ! Look for the extremum of the propagation discriminant
           ! gamma = [g_4*omega**4 + g_2*omega**2 + g_0]/omega**2

           g_4 = -4._WP*V_g*c_1
           g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*l*(l+1)
           g_0 = -4._WP*l*(l+1)*As/c_1

           if(g_0 /= 0._WP .AND. SIGN(1._WP, g_4) == SIGN(1._WP, g_0)) then

              ! Gamma has an extremem (maximum) at omega**2 == omega_ext

              omega2_ext = SQRT(g_0/g_4)
              gamma_ext = (g_4*omega2_ext**2 + g_2*omega2_ext + g_0)/omega2_ext

           else

              ! Gamma has no extremum

              omega2_ext = 0._WP

           endif

           ! Calculate gamma at the a/b frequencies

           gamma_a = (g_4*omega_a**4 + g_2*omega_a**2 + g_0)/omega_a**2
           gamma_b = (g_4*omega_b**4 + g_2*omega_b**2 + g_0)/omega_b**2

           ! Determine the wavenumber bounds

           lambda_osc(i) = SQRT(MAX(-gamma_a, -gamma_b, 0._WP))/x

           if(omega2_ext >= omega_a**2 .AND. omega2_ext <= omega_b**2) then
              lambda_exp(i) = 0.5_WP*(ABS(As + V_g - U + 2._WP - 2._WP*l) + &
                                      SQRT(MAX(gamma_ext, gamma_a, gamma_b, 0._WP)))/x
           else
              lambda_exp(i) = 0.5_WP*(ABS(As + V_g - U + 2._WP - 2._WP*l) + &
                                      SQRT(MAX(gamma_a, gamma_b, 0._WP)))/x
           endif

         end associate
       end associate

    end do wave_loop

    lambda_osc(n) = 0._WP
    lambda_exp(n) = 0._WP

    ! Place points to ensure a given sampling of
    ! oscillatory/exponential scale lengths

    samp_loop : do i = 1,n-1

       dphi_osc = MAX(lambda_osc(i), lambda_osc(i+1))*(x_bc(i+1) - x_bc(i))
       dphi_exp = MAX(lambda_exp(i), lambda_exp(i+1))*(x_bc(i+1) - x_bc(i))

       dn(i) = MAX(dn(i), FLOOR((alpha_osc*dphi_osc)/PI), FLOOR((alpha_exp*dphi_exp)/PI))

    end do samp_loop

    ! Place points to ensure the central evanescent zone has at least
    ! n_center points in it

    ! First, locate the innermost turning point at both omega_a and omega_b

    call find_x_turn(x_bc, bc, op, omega_a, x_turn_a)
    call find_x_turn(x_bc, bc, op, omega_b, x_turn_b)

    x_turn = MIN(x_turn_a, x_turn_b)
    call locate(x_bc, x_turn, i_turn)

    ! Add points

    if(i_turn >= 1 .AND. i_turn < n) then
       dn(:i_turn) = MAX(dn(:i_turn), CEILING(n_center*x_bc(i_turn+1)/x_turn/i_turn))
    endif
    
    ! Place points based on a simple floor
    
    dn = MAX(dn, n_floor)

    ! Finish

    return

  end subroutine plan_dispersion_grid

!****

  subroutine find_x_turn (x, bc, op, omega, x_turn)

    real(WP), intent(in)                     :: x(:)
    class(base_coeffs_t), target, intent(in) :: bc
    type(oscpar_t), target, intent(in)       :: op
    real(WP), intent(in)                     :: omega
    real(WP)                                 :: x_turn

    type(gamma_t) :: gamma
    integer       :: i

    ! Find the inner turning point at frequency omega

    gamma%bc => bc
    gamma%op => op

    gamma%omega = omega

    x_turn = HUGE(0._WP)

    turn_loop : do i = 1,SIZE(x)-1
       if(gamma%eval(x(i)) > 0._WP .AND. gamma%eval(x(i+1)) <= 0._WP) then
          x_turn = gamma%root(x(i), x(i+1), 0._WP)
          exit turn_loop
       end if
    end do turn_loop

    ! Finish

    return

  end subroutine find_x_turn

!****

  function eval_gamma (this, z) result (gamma)

    class(gamma_t), intent(inout) :: this
    complex(WP), intent(in)      :: z
    complex(WP)                  :: gamma

    real(WP) :: x
    real(WP) :: g_4
    real(WP) :: g_2
    real(WP) :: g_0

    ! Calculate the propagation discriminant

    x = REAL(z)

    associate(V_g => this%bc%V(x)/this%bc%Gamma_1(x), As => this%bc%As(x), &
              U => this%bc%U(x), c_1 => this%bc%c_1(x), &
              l => this%op%l)

      g_4 = -4._WP*V_g*c_1
      g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*l*(l+1)
      g_0 = -4._WP*l*(l+1)*As/c_1

      gamma = (g_4*this%omega**4 + g_2*this%omega**2 + g_0)/this%omega**2

    end associate

    ! Finish

    return

  end function eval_gamma

end module gyre_grid
