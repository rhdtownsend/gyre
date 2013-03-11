! Module   : gyre_grid
! Purpose  : grid construction

$include 'core.inc'

module gyre_grid

  ! Uses

  use core_kinds
  use core_constants

  use gyre_mech_coeffs

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

       dx_1 = 0.5_WP/(SUM([(s**(REAL(k-1, WP)/REAL(m-1, WP)),k=1,m)]) - &
                      0.5_WP*s)

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

       dx_1 = 0.5_WP/(SUM([(s**(REAL(k-1, WP)/REAL(m-1, WP)),k=1,m)]))

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

  subroutine plan_dispersion_grid (x_mc, mc, omega, l, alpha_osc, alpha_exp, n_center, n_floor, dn)

    real(WP), intent(in)             :: x_mc(:)
    class(mech_coeffs_t), intent(in) :: mc
    complex(WP), intent(in)          :: omega
    integer, intent(in)              :: l
    real(WP), intent(in)             :: alpha_osc
    real(WP), intent(in)             :: alpha_exp
    integer, intent(in)              :: n_center
    integer, intent(in)              :: n_floor
    integer, intent(inout)           :: dn(:)

    integer     :: i
    complex(WP) :: k_r
    real(WP)    :: dphi_osc
    real(WP)    :: dphi_exp
    real(WP)    :: x_turn

    $CHECK_BOUNDS(SIZE(dn),SIZE(x_mc)-1)

    ! Plan an oversamp grid, modifying dn (see build_grid_oversamp)
    ! with additional points based on a local dispersion analysis

    ! Place points based on the oscillatory (real) and exponential
    ! (imaginary) parts of the local radial wavenumber

    !$OMP PARALLEL DO PRIVATE (k_r, dphi_osc, dphi_exp)
    cell_loop : do i = 1,SIZE(x_mc)-1

       ! Estimate the local radial wavenumber at the cell center

       associate(x_mid => 0.5_WP*(x_mc(i)+x_mc(i+1)))
          associate(V_g => mc%V(x_mid)/mc%Gamma_1(x_mid), As => mc%As(x_mid), c_1 => mc%c_1(x_mid))
            k_r = SQRT(-(l*(l+1)/(c_1*omega**2) - V_g)*(c_1*omega**2 - As))/x_mid
          end associate
       end associate

       ! Place points to ensure a given ratio between cell cell size and
       ! wavelength
       
       dphi_osc = ABS(REAL(k_r))*(x_mc(i+1) - x_mc(i))
       dphi_exp = ABS(AIMAG(k_r))*(x_mc(i+1) - x_mc(i))

       dn(i) = MAX(dn(i), FLOOR((alpha_osc*dphi_osc + alpha_exp*dphi_exp)/PI))

    end do cell_loop

    ! Place points to ensure the central evanescent zone has at least
    ! n_center points in it

    if(l > 0) then
       
       associate(x_2 => x_mc(2))

         associate(V_g => mc%V(x_2)/mc%Gamma_1(x_2), As => mc%As(x_2), c_1 => mc%c_1(x_2))

           if(As > 0._WP) then
              x_turn = MIN(ABS(REAL(SQRT(c_1*omega**2/As))), &
                           ABS(REAL(SQRT(l*(l+1)/(c_1*omega**2*V_g)))))*x_2
           else
              x_turn = ABS(REAL(SQRT(l*(l+1)/(c_1*omega**2*V_g))))*x_2
           endif

         end associate

         dn(1) = MAX(dn(1), CEILING(x_2/x_turn*n_center))

       end associate

    endif

    ! Place points based on a simple floor

    dn = MAX(dn, n_floor)

    ! Finish

    return

  end subroutine plan_dispersion_grid

end module gyre_grid
