! Module   : gyre_grid_weights
! Purpose  : grid weighting schemes and factory procedure
!
! Copyright 2016-2017 Rich Townsend
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

module gyre_grid_weights

  ! Uses

  use core_kinds
  use core_func

  use gyre_grid

  ! No implicit typing

  implicit none

  ! Derived-type definitions (used internally for root-finding)

  type, extends (func_t) :: geom_func_t
     real(WP) :: s
     integer  :: n
   contains
     procedure :: eval_c_ => eval_geom_func_
  end type geom_func_t

  ! Interfaces

  interface grid_t
     module procedure grid_t_weights_
  end interface grid_t

  ! Access specifiers

  private

  public :: grid_t
  public :: uni_weights
  public :: geo_weights
  public :: log_weights

  ! Procedures

contains

  function grid_t_weights_ (w, x_i, x_o) result (gr)

    real(WP), intent(in) :: w(:)
    real(WP), intent(in) :: x_i
    real(WP), intent(in) :: x_o
    type(grid_t)         :: gr

    ! Construct the grid_t using the supplied weights array and range
    ! (x_i,x_o)

    gr = grid_t((1._WP-w)*x_i + w*x_o)

    ! Finish

    return

  end function grid_t_weights_

  !****

  function uni_weights (n) result (w)

    integer, intent(in) :: n
    real(WP)            :: w(n)

    integer :: k

    ! Create an n-point array of weights with uniform spacing across
    ! the [0,1] interval

    weight_loop : do k = 1, n
       w(k) = (k-1._WP)/(n-1._WP)
    end do weight_loop

    ! Finish

    return

  end function uni_weights

  !****

  function geo_weights (n, s) result (w)

    integer, intent(in)   :: n
    real(WP), intent(in)  :: s
    real(WP)              :: w(n)

    integer           :: m
    real(WP)          :: g_a
    real(WP)          :: g_b
    real(WP)          :: g
    type(geom_func_t) :: gf
    real(WP)          :: dw
    integer           :: k

    ! Create an n-point array of weights with geometric spacing in
    ! each half of the [0,1] interval. The parameter s controls
    ! the ratio between the boundary cell size and the average cell
    ! size 1/(n-1)

    if (MOD(n, 2) == 0) then

       ! Even number of points / odd number of cells

       ! Solve for the growth factor g. The upper bound is derived
       ! by applying a Taylor expansion to the equation for g

       m = n/2-1

       gf%n = n
       gf%s = s

       g_a = EPSILON(0._WP)
       g_b = (s*(n-1)-2*m-1)/m

       g = gf%root(g_a, g_b, 0._WP)

       ! Set up weights for the inner part of the interval

       w(1) = 0._WP
       dw = 1._WP/(s*(n-1))
       
       even_weight_loop : do k = 1, m
          w(k+1) = w(k) + dw
          dw = (1._WP+g)*dw
       end do even_weight_loop

       ! Reflect to get the outer part of the interval

       w(m+2:) = 1._WP - w(m+1:1:-1)

    else

       ! Odd number of points / even number of cells

       ! Solve for the growth factor g. The upper bound is derived
       ! by applying a Taylor expansion to the equation for g

       m = (n-1)/2

       gf%n = n
       gf%s = s

       g_a = EPSILON(0._WP)
       g_b = (s*(n-1)-2*m)/(m*(m-1))

       g = gf%root(g_a, g_b, 0._WP)

       ! Set up the inner part of the interval

       w(1) = 0._WP
       dw = 1._WP/(s*(n-1))
       
       odd_weight_loop : do k = 1, m-1
          w(k+1) = w(k) + dw
          dw = (1._WP+g)*dw
       end do odd_weight_loop

       ! Reflect to get the outer part of the interval

       w(m+1:) = 0.5_WP

       w(m+2:) = 1._WP - w(m:1:-1)

    end if

    ! Finish

    return

  end function geo_weights

  !****

  function eval_geom_func_ (this, z) result (f_z)

    class(geom_func_t), intent(inout) :: this
     complex(WP), intent(in)          :: z
     complex(WP)                      :: f_z

     real(WP) :: g
     integer  :: m

     ! Calcuate the discriminant for the geometric growth factor

     g = REAL(z)

     if (MOD(this%n, 2) == 0) then

        m = this%n/2-1

        if (1._WP+g > HUGE(0._WP)**(1._WP/m)) then
           f_z = - (2._WP + g)
        else
           f_z = (2._WP + this%s*(this%n-1)*g)/(1._WP + g)**m - (2._WP + g)
        endif

     else

        m = (this%n-1)/2

        if(1._WP+g > HUGE(0._WP)**(1._WP/m)) then
           f_z = -2._WP
        else
           f_z = (2._WP + this%s*(this%n-1)*g)/(1._WP + g)**m - 2._WP
        endif
        
     endif
     
     ! Finish

    return
    
  end function eval_geom_func_

  !****

  function log_weights (n, s) result (w)

    integer, intent(in)  :: n
    real(WP), intent(in) :: s
    real(WP)             :: w(n)

    real(WP) :: dw_1
    integer  :: k
    real(WP) :: v
    real(WP) :: t

    ! Create an n-point array of weights with geometric spacing in
    ! each half of the [0,1] interval. The parameter s controls
    ! the ratio between the boundary cell size and the average cell
    ! size 1/(n-1)

    dw_1 = 1._WP/(s*(n-1))

    if (MOD(n, 2) == 0) then

       ! Even number of points / odd number of cells

       ! Set up the inner part of the interval

       w(1) = 0._WP

       even_weight_loop : do k = 2, n/2

          v = (k-1.5_WP)/(n/2-1.5_WP)
          t = (1._WP-v)*LOG(0.5_WP) + v*LOG(dw_1)

          w(n/2-k+2) = EXP(t)

       enddo even_weight_loop
       
       ! Reflect to get the outer part of the interval

       w(n/2+1:) = 1._WP - w(n/2:1:-1)

    else

       ! Odd number of points / even number of cells

       ! Set up the inner part of the interval

       w(1) = 0._WP

       odd_weight_loop : do k = 2, (n-1)/2

          v = (k-1._WP)/((n-1)/2-1._WP)
          t = (1._WP-v)*LOG(0.5_WP) + v*LOG(dw_1)

          w((n-1)/2-k+2) = EXP(t)

       end do odd_weight_loop

       w((n+1)/2) = 0.5_WP

       ! Reflect to get the outer part of the interval

       w((n+1)/2+1:) = 1._WP - w((n-1)/2:1:-1)

    end if

    ! Finish

    return

  end function log_weights

end module gyre_grid_weights
