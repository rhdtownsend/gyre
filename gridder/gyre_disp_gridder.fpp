! Module   : gyre_disp_gridder
! Purpose  : dispersion analysis grid construction

$include 'core.inc'

module gyre_disp_gridder

  ! Uses

  use core_kinds

  use gyre_gridder
  use gyre_mech_coeffs

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(gridder_t) :: disp_gridder_t
     private
     class(mech_coeffs_t), pointer :: mc => null()
     real(WP), allocatable         :: x(:)
     real(WP)                      :: f_wave_osc
     real(WP)                      :: f_wave_exp
     integer                       :: n_center
     integer                       :: n_floor
     integer                       :: l
     integer                       :: n
   contains
     private
     procedure, public :: init
     procedure, public :: build
  end type disp_gridder_t

  ! Access specifiers

  private

  public :: disp_gridder_t

  ! Procedures

contains

  subroutine init (this, x, mc, f_wave_osc, f_wave_exp, n_center, n_floor, l)

    class(disp_gridder_t), intent(out)       :: this
    real(WP), intent(in)                     :: x(:)
    class(mech_coeffs_t), intent(in), target :: mc
    real(WP), intent(in)                     :: f_wave_osc
    real(WP), intent(in)                     :: f_wave_exp
    integer, intent(in)                      :: n_center
    integer, intent(in)                      :: n_floor
    integer, intent(in)                      :: l
   
    ! Initialize the disp_gridder

    this%x = x
    this%mc => mc
    
    this%f_wave_osc = f_wave_osc
    this%f_wave_evan = f_wave_evan
    this%n_center = n_center
    this%n_floor = n_floor

    this%l = l
    this%n = SIZE(x)

    ! Finish

    return

  end subroutine init

!****

  subroutine build (this, omega, x)

    class(disp_gridder_t), intent(in)  :: this
    complex(WP), intent(in)            :: omega
    real(WP), allocatable, intent(out) :: x(:)

    integer :: n_add(this%n-1)
    integer :: n
    integer :: i
    integer :: j
    integer :: k

    ! Build a grid based on oversampling this%x, with additional
    ! points based on a local dispersion analysis

    ! First, determine the number of additional points to place in
    ! each cell

    n_add = 0

    call this%add_wave_points(omega, n_add)
    call this%add_center_points(omega, n_add)
    call this%add_floor_points(n_add)

    ! Now set up the grid

    n = SIZE(this%x)

    allocate(x(SUM(n_add) + n))

    k = 1

    do i = 1,n-1
       do j = 1,n_add(i)+1
          x(k) = this%x(i) + (j-1)*(this%x(i+1)-this%x(i))/(n_add(i)+1)
          k = k + 1
       end do
    end do
    
    x(k) = this%x(n)

    ! Finish

    return

  end subroutine build

!****

  subroutine add_wave_points (this, omega, n_add)

    class(disp_gridder_t), intent(in) :: this
    complex(WP), intent(in)           :: omega
    integer, intent(inout)            :: n_add(:)

    integer     :: i
    complex(WP) :: k_r
    real(WP)    :: dphi_osc
    real(WP)    :: dphi_exp

    $CHECK_BOUNDS(SIZE(n_add),this%n-1)

    ! Add points based on the oscillatory (real) and exponential
    ! (imaginary) parts of the local radial wavenumber

    cell_loop : do i = 1,this%n-1

       ! Estimate the local radial wavenumber at the cell center

       associate(x_mid => 0.5_WP*(x_cf(i)+x_cf(i+1)))
          associate(V_g => this%cf%V(x_mid)/this%cf%Gamma_1(x_mid), &
                    As => this%cf%As(x_mid), c_1 => this%cf%c_1(x_mid), &
                    l => this%l)
            k_r = SQRT(-(l*(l+1)/(c_1*omega**2) - V_g)*(c_1*omega**2 - As))/x_mid
          end associate
       end associate

       ! Add points

       dphi_osc = ABS(REAL(k_r))*(cf%x(i+1) - cf%x(i))
       dphi_exp = ABS(AIMAG(k_r))*(cf%x(i+1) - cf%x(i))

       n_add(i) = MAX(n_add(i), FLOOR((this%f_wave_osc*dphi_osc + this%f_wave_exp*dphi_exp)/PI))

    end do count_loop

    ! Finish

    return

  end subroutine add_wave_points

!****

  subroutine add_center_points (this, omega, n_add)

    class(disp_gridder_t), intent(in) :: this
    complex(WP), intent(in)           :: omega
    integer, intent(inout)            :: n_add(:)

    $CHECK_BOUNDS(SIZE(n_add),this%n-1)

    ! Add points to ensure the central evanescent zone has at least
    ! n_center points in it

    if(this%l > 0) then
       
       associate(x_2 => this%x(2))

         associate(V_g => cf%V(x_2)/ac%Gamma_1(x_2), As => cf%As(x_2), c_1 => cf%c_1(x_2))

           if(As > 0._WP) then
              x_turn = SQRT(MIN(c_1*omega_ref**2/As, l_ref*(l_ref+1)/(c_1*omega_ref**2*V_g)))*x_2
           else
              x_turn = SQRT(l_ref*(l_ref+1)/(c_1*omega_ref**2*V_g))*x_2
           endif

         end associate

         n_add(1) = MAX(n_add(1), CEILING(x_2/x_turn*this%n_center))

       end associate

    endif

    ! Finish

    return

  end subroutine add_center_points

!****

  subroutine add_floor_points (this, n_add)

    class(disp_gridder_t), intent(in) :: this
    integer, intent(inout)            :: n_add(:)

    $CHECK_BOUNDS(SIZE(n_add),this%n-1)

    ! Add points based on a simple floor

    n_add = MAX(n_add, this%n_floor)

    ! Finish

    return

  end subroutine add_floor_points

end module gyre_disp_gridder
