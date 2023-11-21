! Module  : lane_emden_m
! Purpose : Lane-Emden equation solver for composite polytropes
!
! Copyright 2015-2020 Rich Townsend & The GYRE Team
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

module lane_emden_m

  ! Uses

  use kinds_m
  use memory_m

  use math_m

  use odepack

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  real(WP), save :: n_poly_m
  real(WP), save :: B_m
  real(WP), save :: x_m

  !$OMP THREADPRIVATE(n_poly_m, B_m, x_m)

  ! Access specifiers

  private

  public :: solve_lane_emden

  ! Procedures

contains

  subroutine solve_lane_emden (n_poly, z_b, Delta_b, dz, tol, z, theta, dtheta, n_r_out)

    real(WP), intent(in)               :: n_poly(:)
    real(WP), intent(in)               :: z_b(:)
    real(WP), intent(in)               :: Delta_b(:)
    real(WP), intent(in)               :: dz
    real(WP), intent(in)               :: tol
    real(WP), allocatable, intent(out) :: z(:)
    real(WP), allocatable, intent(out) :: theta(:)
    real(WP), allocatable, intent(out) :: dtheta(:)
    integer, intent(out)               :: n_r_out

    real(WP), parameter :: X_BEG = sqrt(EPSILON(0._WP))
    integer, parameter  :: D_0 = 512
    integer, parameter  :: NEQ = 2
    integer, parameter  :: NG = 2

    integer               :: d
    integer               :: n
    real(WP), allocatable :: x(:)
    real(WP), allocatable :: y(:,:)
    real(WP)              :: dx
    real(WP)              :: x_exp
    real(WP)              :: y_exp(2)
    integer               :: n_r
    logical               :: first
    integer               :: i
    integer               :: istate
    real(WP)              :: rwork(22+NEQ*MAX(16,NEQ+9)+3*NG)
    integer               :: iwork(20+NEQ)
    integer               :: jroot(NG)
    real(WP)              :: t_t
    real(WP)              :: f

    $CHECK_BOUNDS(SIZE(z_b),SIZE(n_poly)-1)
    $CHECK_BOUNDS(SIZE(Delta_b),SIZE(n_poly)-1)

    if (SIZE(z_b) > 0) then
       $ASSERT(z_b(1) > X_BEG,First boundary too close to origin)
    endif

    ! Initialize arrays

    d = D_0

    allocate(x(d))
    allocate(y(2,d))

    n = 0

    ! Perform a series expansion out to X_BEG

    dx = dz

    n_poly_m = n_poly(1)
    B_m = 1._WP

    x_exp = 0._WP

    expand_loop : do

       ! Use a forth-order expansion of the Lane-Emden solutions

       y_exp(1) = 1._WP - x_exp**2/6._WP + n_poly_m*x_exp**4/120._WP
       y_exp(2) = -x_exp/3._WP + n_poly_m*x_exp**3/30._WP

       if (x_exp < X_BEG) then

          ! If necessary, expand arrays

          n = n + 1

          if(n > d) then
             d = 2*d
             call reallocate(x, [d])
             call reallocate(y, [2,d])
          endif

          x(n) = x_exp
          y(:,n) = y_exp

       else

          exit expand_loop

       end if

       x_exp = MIN(x_exp+dx, X_BEG)

    end do expand_loop

    ! Now continue by intergrating to each boundary point, or to the
    ! surface

    n_r = SIZE(z_b) + 1

    first = .TRUE.

    region_loop : do i = 1, n_r

       n_poly_m = n_poly(i)

       if (i < n_r) then
          x_m = z_b(i)
       else
          x_m = HUGE(0._WP)
       endif

       istate = 1

       integrate_loop : do

          ! Integrate one step

          n = n + 1

          if(n > d) then
             d = 2*d
             call reallocate(x, [d])
             call reallocate(y, [2,d])
          endif

          if (first) then
             x(n) = x_exp
             y(:,n) = y_exp
             first = .FALSE.
          else
             x(n) = x(n-1)
             y(:,n) = y(:,n-1)
          end if

          call LSODAR(lane_emden_rhs, [NEQ], y(:,n), x(n), x(n-1)+dx, 1, [0._WP,0._WP], [tol,tol], &
               1, istate, 0, rwork, SIZE(rwork), iwork, SIZE(iwork), lane_emden_jac, 1, &
               lane_emden_constr, NG, jroot)

          ! Check for success

          select case (istate)
          case (2)
          case (3)
             if (jroot(1) == 1) then
                exit integrate_loop
             elseif (jroot(2) == 1) then
                exit region_loop
             else
                $ABORT(Root not found)
             endif
          case default
             write(ERROR_UNIT, *) 'istate is ',istate
             write(ERROR_UNIT, *) 'Integration range x=',x(n-1)
             write(ERROR_UNIT, *) 'Starting location y=',y(:,n-1)
             write(ERROR_UNIT, *) 'Relative tolerance =',0._WP
             write(ERROR_UNIT, *) 'Absolute tolerance =',tol
             write(ERROR_UNIT, *) 'Number of steps    =',n
             $ABORT(Integration failed)
          end select
       
       end do integrate_loop

       ! Create the boundary double point

       n = n + 1

       if(n > d) then
          d = 2*d
          call reallocate(x, [d])
          call reallocate(y, [2,d])
       endif

       x(n) = x(n-1)

       t_t = exp(n_poly(i)*log(y(1,n-1)) + Delta_b(i))

       y(1,n) = 1._WP

       f = ((n_poly(i  )+1._WP)*y(1,n  )**(n_poly(i+1)+1._WP))/ &
           ((n_poly(i+1)+1._WP)*y(1,n-1)**(n_poly(i  )+1._WP))

       y(2,n) = f*t_t*y(2,n-1)

       ! Update B 

       B_m = f*t_t**2*B_m

    end do region_loop

    y(1,n) = 0._WP

    ! Record how many regions were actually encountered

    n_r_out = i

    ! Set up return values

    z = x(:n)
    
    theta = y(1,:n)
    dtheta = y(2,:n)

    ! Finish

    return

  end subroutine solve_lane_emden

!****

  subroutine lane_emden_rhs (neq, x, y, dy_dx)

    integer, intent(in)   :: neq
    real(WP), intent(in)  :: x
    real(WP), intent(in)  :: y(neq)
    real(WP), intent(out) :: dy_dx(neq)

    ! Calculate the right-hand side vector

    dy_dx(1) = y(2)
    if (n_poly_m /= 0._WP) then
       dy_dx(2) = -B_m*abs(y(1))**n_poly_m - 2._WP*y(2)/x
    else
       dy_dx(2) = -B_m - 2._WP*y(2)/x
    endif

    ! Finish

    return

  end subroutine lane_emden_rhs

!****

  subroutine lane_emden_jac (neq, x, y, ml, mu, pd, nrowpd)

    integer, intent(in)   :: neq
    real(WP), intent(in)  :: x
    real(WP), intent(in)  :: y(neq)
    integer, intent(in)   :: ml
    integer, intent(in)   :: mu
    integer, intent(in)   :: nrowpd
    real(WP), intent(out) :: pd(nrowpd,neq)

    ! Calculate the Jacobian matrix

    pd(1,1) = 0._WP
    pd(1,2) = 1._WP

    pd(2,1) = -n_poly_m*B_m*abs(y(1))**(n_poly_m-1._WP)*SIGN(1._WP, y(1))
    pd(2,2) = -2._WP/x

    ! Finish

    return

  end subroutine lane_emden_jac

!****

  subroutine lane_emden_constr (neq, x, y, ng, gout)

    integer, intent(in)   :: neq
    real(WP), intent(in)  :: x
    real(WP), intent(in)  :: y(neq)
    integer, intent(in)   :: ng
    real(WP), intent(out) :: gout(ng)

    ! Calculate the constraint function

    gout(1) = x - x_m ! Region boundary
    gout(2) = y(1)    ! Surface

    ! Finish

    return

  end subroutine lane_emden_constr

end module lane_emden_m
