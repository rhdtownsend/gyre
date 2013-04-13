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
  use gyre_gridpar

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

  public :: build_grid

  ! Procedures

contains

  subroutine build_grid (gp, bc, op, x_in, x)

    type(gridpar_t), intent(in)        :: gp(:)
    class(base_coeffs_t), intent(in)   :: bc
    type(oscpar_t), intent(in)         :: op
    real(WP), allocatable, intent(in)  :: x_in(:)
    real(WP), allocatable, intent(out) :: x(:)

    integer :: i

    ! Build a grid using the supplied list of gridpars

    select case (gp(1)%op_type)
    case ('CREATE_CLONE')
       x = x_in
    case ('CREATE_GEOM')
       call create_geom(gp(1)%s, gp(1)%n, x)
    case ('CREATE_LOG')
       call create_log(gp(1)%s, gp(1)%n, x)
    case default
       $ABORT(Invalid op_type (the first op_type must be CREATE_*))
    end select

    gp_loop : do i = 2,SIZE(gp)

       select case (gp(i)%op_type)
       case ('RESAMP_DISPERSION')
          call resample_dispersion(bc, op, gp(i)%omega_a, gp(i)%omega_b, &
                                   gp(i)%alpha_osc, gp(i)%alpha_exp, x)
       case ('RESAMP_CENTER')
          call resample_center(bc, op, gp(i)%omega_a, gp(i)%omega_b, gp(i)%n, x)
       case ('RESAMP_UNIFORM')
          call resample_uniform(gp(i)%n, x)
       case default
          $ABORT(Invalid op_type (the second and subsequent op_types must be RESAMPLE_*))
       end select

    end do gp_loop

    ! Finish

    return

  end subroutine build_grid
          
!****

  subroutine create_geom (s, n, x)

    real(WP), intent(in)               :: s
    integer, intent(in)                :: n
    real(WP), allocatable, intent(out) :: x(:)

    integer  :: m
    integer  :: k
    real(WP) :: dx_1

    ! Create an n-point grid with geometric spacing in each half of the
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

  end subroutine create_geom

!****

  subroutine create_log (s, n, x)

    real(WP), intent(in)               :: s
    integer, intent(in)                :: n
    real(WP), allocatable, intent(out) :: x(:)

    real(WP) :: dx_1
    integer  :: k
    real(WP) :: w
    real(WP) :: t

    ! Create an n-point grid with logarithmic spacing in each half of
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

  end subroutine create_log

!****

  subroutine resample (x, dn)

    real(WP), allocatable, intent(inout) :: x(:)
    integer, intent(in)                  :: dn(:)

    integer               :: n
    real(WP), allocatable :: x_new(:)
    integer               :: i
    integer               :: j
    integer               :: k
    
    $CHECK_BOUNDS(SIZE(dn),SIZE(x)-1)

    ! Resample x by inserting dn additional points across the cells
    
    n = SIZE(x)

    allocate(x_new(SUM(dn) + n))

    k = 1

    do i = 1,n-1
       do j = 1,dn(i)+1
          x_new(k) = x(i) + (j-1)*(x(i+1)-x(i))/(dn(i)+1)
          k = k + 1
       end do
    end do
    
    x_new(k) = x(n)

    call MOVE_ALLOC(x_new, x)

    ! Finish

    return

  end subroutine resample

!****

  subroutine resample_dispersion (bc, op, omega_a, omega_b, alpha_osc, alpha_exp, x)

    class(base_coeffs_t), intent(in)     :: bc
    type(oscpar_t), intent(in)           :: op
    real(WP), intent(in)                 :: omega_a
    real(WP), intent(in)                 :: omega_b
    real(WP), intent(in)                 :: alpha_osc
    real(WP), intent(in)                 :: alpha_exp
    real(WP), allocatable, intent(inout) :: x(:)

    integer               :: n_x
    real(WP), allocatable :: k_osc(:)
    real(WP), allocatable :: k_exp(:)
    integer               :: i
    real(WP)              :: g_4
    real(WP)              :: g_2
    real(WP)              :: g_0
    real(WP)              :: omega2_ext
    real(WP)              :: gamma_ext
    real(WP)              :: gamma_a
    real(WP)              :: gamma_b
    integer, allocatable  :: dn(:)
    real(WP)              :: dphi_osc
    real(WP)              :: dphi_exp

    $ASSERT(ALLOCATED(x),No input grid)

    $ASSERT(omega_b >= omega_a,Invalid frequency interval)

    ! Resample x by adding points to each cell, such that there are at
    ! least alpha_osc points per oscillatory wavelength and alpha_exp
    ! points per exponential wavelength. Wavelengths are calculated
    ! based on a local dispersion analysis of the adibatic/Cowling
    ! wave equation, for frequencies in the interval [omega_a,omega_b]

    ! First, determine maximal oscillatory and exponential wavenumbers
    ! at each point of x

    n_x = SIZE(x)

    allocate(k_osc(n_x))
    allocate(k_exp(n_x))

    k_osc(1) = 0._WP
    k_exp(1) = 0._WP

    wavenumber_loop : do i = 2,n_x-1

       associate(V_g => bc%V(x(i))/bc%Gamma_1(x(i)), As => bc%As(x(i)), &
                 U => bc%U(x(i)), c_1 => bc%c_1(x(i)), &
                 l => op%l)

         ! Look for an extremum of the propagation discriminant ]
         ! gamma = [g_4*omega**4 + g_2*omega**2 + g_0]/omega**2 in the
         ! frequency interval

         g_4 = -4._WP*V_g*c_1
         g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*l*(l+1)
         g_0 = -4._WP*l*(l+1)*As/c_1

         if(g_0 /= 0._WP .AND. SIGN(1._WP, g_4) == SIGN(1._WP, g_0)) then

            omega2_ext = SQRT(g_0/g_4)
            
            if(omega2_ext >= omega_a**2 .AND. omega2_ext <= omega_b**2) then
               gamma_ext = (g_4*omega2_ext**2 + g_2*omega2_ext + g_0)/omega2_ext
            else
               gamma_ext = 0._WP
            endif

         else

            gamma_ext = 0._WP

         endif

         ! Calculate gamma at the interval endpoints

         gamma_a = (g_4*omega_a**4 + g_2*omega_a**2 + g_0)/omega_a**2
         gamma_b = (g_4*omega_b**4 + g_2*omega_b**2 + g_0)/omega_b**2

         ! Determine the wavenumber bounds

         k_osc(i) = SQRT(MAX(-gamma_a, -gamma_b, 0._WP))/x(i)

         k_exp(i) = 0.5_WP*(ABS(As + V_g - U + 2._WP - 2._WP*l) + &
                            SQRT(MAX(gamma_ext, gamma_a, gamma_b, 0._WP)))/x(i)

       end associate

    end do wavenumber_loop

    k_osc(n_x) = 0._WP
    k_exp(n_x) = 0._WP

    ! Determine how many points to add to each cell

    allocate(dn(n_x-1))

    sample_loop : do i = 1,n_x-1

       ! Calculate the oscillatory and exponential phase change across
       ! the cell

       dphi_osc = MAX(k_osc(i), k_osc(i+1))*(x(i+1) - x(i))
       dphi_exp = MAX(k_exp(i), k_exp(i+1))*(x(i+1) - x(i))

       ! Set up dn

       dn(i) = MAX(FLOOR((alpha_osc*dphi_osc)/TWOPI), FLOOR((alpha_exp*dphi_exp)/TWOPI))

    end do sample_loop

    ! Perform the resampling

    call resample(x, dn)

    ! Finish

    return

  end subroutine resample_dispersion

!****

  subroutine resample_center (bc, op, omega_a, omega_b, n, x)

    class(base_coeffs_t), intent(in)     :: bc
    type(oscpar_t), intent(in)           :: op
    real(WP), intent(in)                 :: omega_a
    real(WP), intent(in)                 :: omega_b
    integer, intent(in)                  :: n
    real(WP), allocatable, intent(inout) :: x(:)

    real(WP)             :: x_turn_a
    real(WP)             :: x_turn_b
    real(WP)             :: x_turn
    integer              :: i_turn
    integer              :: n_x
    integer, allocatable :: dn(:)

    $ASSERT(ALLOCATED(x),No input grid)

    $ASSERT(omega_b >= omega_a,Invalid frequency interval)

    ! Resample x by adding points to central cells, such that there
    ! are n_floor points covering the evanescent region at the
    ! center. Evanescence is determined based on a local dispersion
    ! analysis of the adibatic/Cowling wave equation, for frequencies
    ! in the interval [omega_a,omega_b]

    ! First, locate the innermost turning point at both omega_a and omega_b

    call find_x_turn(x, bc, op, omega_a, x_turn_a)
    call find_x_turn(x, bc, op, omega_b, x_turn_b)

    x_turn = MIN(x_turn_a, x_turn_b)
    call locate(x, x_turn, i_turn)

    ! Determine how many points to add to each cell

    n_x = SIZE(x)

    allocate(dn(n_x-1))

    dn = 0

    if(i_turn >= 1 .AND. i_turn < n_x) then
       dn(:i_turn) = CEILING(n*x(i_turn+1)/x_turn/i_turn)
    endif
    
    ! Perform the resampling

    call resample(x, dn)

    ! Finish

    return

  end subroutine resample_center

!****

  subroutine resample_uniform (n, x)

    integer, intent(in)                  :: n
    real(WP), allocatable, intent(inout) :: x(:)

    integer              :: n_x
    integer, allocatable :: dn(:)

    $ASSERT(ALLOCATED(x),No input grid)

    ! Resample x by adding n points to each cell

    n_x = SIZE(x)

    allocate(dn(n_x-1))

    dn = n

    ! Perform the resampling

    call resample(x, dn)

    ! Finish

    return

  end subroutine resample_uniform

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
