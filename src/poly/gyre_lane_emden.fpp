! Module   : gyre_lane_emden
! Purpose  : Lane-Emden equation solver with discontinuities
!
! Copyright 2015-2016 Rich Townsend
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

module gyre_lane_emden

  ! Uses

  use core_kinds
  use core_memory

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

  subroutine solve_lane_emden (n_poly, xi_d, Delta_d, dxi, tol, xi, Theta, dTheta)

    real(WP), intent(in)               :: n_poly(:)
    real(WP), intent(in)               :: xi_d(:)
    real(WP), intent(in)               :: Delta_d(:)
    real(WP), intent(in)               :: dxi
    real(WP), intent(in)               :: tol
    real(WP), allocatable, intent(out) :: xi(:)
    real(WP), allocatable, intent(out) :: Theta(:)
    real(WP), allocatable, intent(out) :: dTheta(:)

    real(WP), parameter :: X_BEG = SQRT(EPSILON(0._WP))
    integer, parameter  :: D_0 = 512
    integer, parameter  :: NEQ = 2
    integer, parameter  :: NG = 1

    integer               :: d
    integer               :: n
    real(WP), allocatable :: x(:)
    real(WP), allocatable :: y(:,:)
    real(WP)              :: dx
    integer               :: n_d
    integer               :: i
    integer               :: istate
    real(WP)              :: rwork(22+NEQ*MAX(16,NEQ+9)+3*NG)
    integer               :: iwork(20+NEQ)
    integer               :: jroot(NG)
    real(WP)              :: t_t
    real(WP)              :: f

    $CHECK_BOUNDS(SIZE(xi_d),SIZE(n_poly)-1)
    $CHECK_BOUNDS(SIZE(Delta_d),SIZE(n_poly)-1)

    ! Initialize arrays

    d = D_0

    allocate(x(d))
    allocate(y(2,d))

    n = 0

    ! Perform a series expansion out to x_beg

    dx = dxi

    n_poly_m = n_poly(1)
    B_m = 1._WP

    expand_loop : do

       ! If necessary, expand arrays

       n = n + 1

       if(n > d) then
          d = 2*d
          call reallocate(x, [d])
          call reallocate(y, [2,d])
       endif

       ! Use a forth-order expansion of the Lane-Emden solutions

       x(n) = dx*(n-1)

       y(1,n) = 1._WP - x(n)**2/6._WP + n_poly_m*x(n)**4/120._WP
       y(2,n) = -x(n)/3._WP + n_poly_m*x(n)**3/30._WP

       if(x(n) > X_BEG) exit expand_loop

    end do expand_loop

    ! Now continue by intergrating to each discontinuity point

    n_d = SIZE(xi_d)

    do i = 1, n_d

       n_poly_m = n_poly(i)
       x_m = xi_d(i)

       istate = 1

       integrate_d_loop : do

          ! Integrate one step

          n = n + 1

          if(n > d) then
             d = 2*d
             call reallocate(x, [d])
             call reallocate(y, [2,d])
          endif

          x(n) = x(n-1)
          y(:,n) = y(:,n-1)

          call LSODAR(lane_emden_rhs, [NEQ], y(:,n), x(n), x(n-1)+dx, 1, [0._WP,0._WP], [tol,tol], &
               1, istate, 0, rwork, SIZE(rwork), iwork, SIZE(iwork), lane_emden_jac, 1, &
               lane_emden_constr_d, NG, jroot)

          ! Check for success

          select case (istate)
          case (2)
          case (3)
             $ASSERT(jroot(1)==1,Root not found)
             exit integrate_d_loop
          case default
             write(ERROR_UNIT, *) 'istate is ',istate
             write(ERROR_UNIT, *) 'Integration range x=',x(n-1)
             write(ERROR_UNIT, *) 'Starting location y=',y(:,n-1)
             write(ERROR_UNIT, *) 'Relative tolerance =',0._WP
             write(ERROR_UNIT, *) 'Absolute tolerance =',tol
             write(ERROR_UNIT, *) 'Number of steps    =',n
             $ABORT(Integration failed)
          end select
       
       end do integrate_d_loop

       ! Create the discontinuity double point

       n = n + 1

       if(n > d) then
          d = 2*d
          call reallocate(x, [d])
          call reallocate(y, [2,d])
       endif

       x(n) = x(n-1)

       t_t = EXP(n_poly(i)*LOG(y(1,n-1)) + Delta_d(i))

       y(1,n) = 1._WP

       f = ((n_poly(i  )+1._WP)*y(1,n  )**(n_poly(i+1)+1._WP))/ &
           ((n_poly(i+1)+1._WP)*y(1,n-1)**(n_poly(i  )+1._WP))

       y(2,n) = f*t_t*y(2,n-1)

       ! Update B 

       B_m = f*t_t**2*B_m

    end do

    ! Finish by integrating to the surface

    n_poly_m = n_poly(n_d+1)

    istate = 1

    integrate_s_loop : do

       ! Integrate one step

       n = n + 1

       if(n > d) then
          d = 2*d
          call reallocate(x, [d])
          call reallocate(y, [2,d])
       endif
       
       x(n) = x(n-1)
       y(:,n) = y(:,n-1)
       
       call LSODAR(lane_emden_rhs, [NEQ], y(:,n), x(n), x(n-1)+dx, 1, [0._WP,0._WP], [tol,tol], &
            1, istate, 0, rwork, SIZE(rwork), iwork, SIZE(iwork), lane_emden_jac, 1, &
            lane_emden_constr_S, NG, jroot)

       ! Check for success

       select case (istate)
       case (2)
       case (3)
          $ASSERT(jroot(1)==1,Root not found)
          exit integrate_s_loop
       case default
          write(ERROR_UNIT, *) 'istate is ',istate
          write(ERROR_UNIT, *) 'Integration range x=',x(n-1)
          write(ERROR_UNIT, *) 'Starting location y=',y(:,n-1)
          write(ERROR_UNIT, *) 'Relative tolerance =',0._WP
          write(ERROR_UNIT, *) 'Absolute tolerance =',tol
          write(ERROR_UNIT, *) 'Number of steps    =',n
          $ABORT(Integration failed)
       end select
       
    end do integrate_s_loop

    y(1,n) = 0._WP

    ! Set up return values

    xi = x(:n)
    
    Theta = y(1,:n)
    dTheta = y(2,:n)

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
       dy_dx(2) = -B_m*ABS(y(1))**n_poly_m - 2._WP*y(2)/x
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

    ! Calculate the Jacobian matrix (segment I)

    pd(1,1) = 0._WP
    pd(1,2) = 1._WP

    pd(2,1) = -n_poly_m*B_m*ABS(y(1))**(n_poly_m-1._WP)*SIGN(1._WP, y(1))
    pd(2,2) = -2._WP/x

    ! Finish

    return

  end subroutine lane_emden_jac

!****

  subroutine lane_emden_constr_d (neq, x, y, ng, gout)

    integer, intent(in)   :: neq
    real(WP), intent(in)  :: x
    real(WP), intent(in)  :: y(neq)
    integer, intent(in)   :: ng
    real(WP), intent(out) :: gout(ng)

    ! Calculate the constraint function (discontinuity)

    gout(1) = x - x_m

    ! Finish

    return

  end subroutine lane_emden_constr_d

!****

  subroutine lane_emden_constr_S (neq, x, y, ng, gout)

    integer, intent(in)   :: neq
    real(WP), intent(in)  :: x
    real(WP), intent(in)  :: y(neq)
    integer, intent(in)   :: ng
    real(WP), intent(out) :: gout(ng)

    ! Calculate the constraint function (surface)

    gout(1) = y(1)

    ! Finish

    return

  end subroutine lane_emden_constr_S

end module gyre_lane_emden
