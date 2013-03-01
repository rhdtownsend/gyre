! Module   : gyre_log_gridder
! Purpose  : log grid construction

$include 'core.inc'

module gyre_log_gridder

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(gridder_t) :: log_gridder_t
     private
     integer  :: n
     real(WP) :: s
   contains
     private
     procedure, public :: init
     procedure, public :: build
  end type log_gridder_t

  ! Access specifiers

  private

  public :: log_gridder_t

  ! Procedures

contains

  subroutine init (this, n, s)

    class(log_gridder_t), intent(out) :: this
    integer, intent(in)               :: n
    real(WP), intent(in)              :: s

    ! Initialize the log_gridder

    this%n = n
    this%s = s

    ! Finish

    return

  end subroutine init

!****

  subroutine build (this, omega, x)

    class(log_gridder_t), intent(in)   :: this
    complex(WP), intent(in)            :: omega
    real(WP), intent(out), allocatable :: x(:)

    real(WP) :: dx_1
    integer  :: k
    real(WP) :: w
    real(WP) :: t

    ! Build a grid with logarithmic spacing in each half of the
    ! [0,1] interval. The parameter s controls the ratio between the
    ! mean cell size and the boundary cell size. omega is not used

    allocate(x(this%n))

    dx_1 = 1._WP/(this%s*(this%n-1))

    if(MOD(this%n, 2) == 0) then

       ! Even number of grid points / odd number of cells

       ! Set up the inner part of the grid

       x(1) = 0._WP

       even_grid_loop : do k = 2,this%n/2

          w = (k-1.5_WP)/(this%n/2-1.5_WP)
          t = LOG(0.5_WP)*(1._WP-w) + LOG(dx_1)*w

          x(this%n/2-k+2) = EXP(t)

       enddo even_grid_loop
       
       ! Reflect to get the outer part of the grid

       x(this%n/2+1:) = 1._WP - x(this%n/2:1:-1)

    else

       ! Odd number of grid points / even number of cells

       ! Set up the inner part of the grid

       x(1) = 0._WP

       odd_grid_loop : do k = 2,(this%n-1)/2

          w = (k-1._WP)/((this%n-1)/2-1._WP)
          t = LOG(0.5_WP)*(1._WP-w) + LOG(dx_1)*w

          x((this%n-1)/2-k+2) = EXP(t)

       end do odd_grid_loop

       x((this%n+1)/2) = 0.5_WP

       ! Reflect to get the outer part of the grid

       x((this%n+1)/2+1:) = 1._WP - x((this%n-1)/2:1:-1)

    end if

    ! Finish

    return

  end subroutine build

end module gyre_log_gridder
