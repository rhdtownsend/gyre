! Module   : gyre_geom_gridder
! Purpose  : geom grid construction

$include 'core.inc'

module gyre_geom_gridder

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(gridder_t) :: geom_gridder_t
     private
     integer  :: n
     real(WP) :: s
   contains
     private
     procedure, public :: init
     procedure, public :: build
  end type geom_gridder_t

  ! Access specifiers

  private

  public :: geom_gridder_t

  ! Procedures

contains

  subroutine init (this, n, s)

    class(geom_gridder_t), intent(out) :: this
    integer, intent(in)                :: n
    real(WP), intent(in)               :: s

    ! Initialize the geom_gridder

    this%n = n
    this%s = s

    ! Finish

    return

  end subroutine init

!****

  subroutine build (this, omega, x)

    class(geom_gridder_t), intent(in)   :: this
    complex(WP), intent(in)             :: omega
    real(WP), intent(out), allocatable  :: x(:)

    integer  :: m
    integer  :: k
    real(WP) :: dx_1
 
    ! Build a grid with geometric spacing in each half of the [0,1]
    ! interval. The parameter s controls the ratio between the central
    ! cell size and the boundary cell size. omega is not used

    allocate(x(this%n))

    if(MOD(this%n, 2) == 0) then

       ! Even number of grid points / odd number of cells

       ! Solve for the size of the boundary cells

       m = this%n/2

       dx_1 = 0.5_WP/(SUM([(this%s**(REAL(k-1, WP)/REAL(m-1, WP)),k=1,m)]) - &
                      0.5_WP*this%s)

       ! Set up the inner part of the grid

       x(1) = 0._WP
       
       even_grid_loop : do k = 2,m
          x(k) = x(k-1) + dx_1*this%s**(REAL(k-2, WP)/REAL(m-1, WP))
       end do even_grid_loop

       ! Reflect to get the outer part of the grid

       x(m+1:) = 1._WP - x(m:1:-1)

    else

       ! Odd number of grid points / even number of cells

       ! Solve for the size of the boundary cells

       m = (this%n-1)/2

       dx_1 = 0.5_WP/(SUM([(this%s**(REAL(k-1, WP)/REAL(m-1, WP)),k=1,m)]))

       ! Set up the inner part of the grid

       x(1) = 0._WP
       
       odd_grid_loop : do k = 2,m
          x(k) = x(k-1) + dx_1*this%s**(REAL(k-2, WP)/REAL(m-1, WP))
       end do odd_grid_loop

       x(m+1) = 0.5_WP

       ! Reflect to get the outer part of the grid

       x(m+2:) = 1._WP - x(m:1:-1)

    end if

    ! Finish

    return

  end subroutine build

end module gyre_geom_gridder
