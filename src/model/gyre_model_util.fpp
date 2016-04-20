! Module   : gyre_model_util
! Purpose  : stellar model utilities
!
! Copyright 2016 Rich Townsend
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

module gyre_model_util

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_model_par
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: set_uniform_rot
  public :: snap_points

  ! Procedures

contains

  !****

  subroutine set_uniform_rot (ml_p, M_star, R_star, Omega_rot)

    type(model_par_t), intent(in) :: ml_p
    real(WP), intent(in)          :: M_star
    real(WP), intent(in)          :: R_star
    real(WP), intent(out)         :: Omega_rot(:)

    ! Set the uniform dimensionless rotation rate

    select case (ml_p%Omega_units)
    case ('HZ')
       Omega_rot = TWOPI*ml_p%Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    case ('UHZ')
       Omega_rot = TWOPI*ml_p%Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))/1E6
    case ('PER_DAY')
       Omega_rot = TWOPI*ml_p%Omega_rot/86400._WP
    case ('CRITICAL')
       Omega_rot = ml_p%Omega_rot*SQRT(8._WP/27._WP)
    case default
       Omega_rot = ml_p%Omega_rot
    end select

    ! Finish

    return

  end subroutine set_uniform_rot
    
  !****

  subroutine snap_points (dx_snap, x, m)

    real(WP), intent(in)              :: dx_snap
    real(WP), intent(inout)           :: x(:)
    real(WP), optional, intent(inout) :: m(:)

    integer  :: i
    real(WP) :: x_snap
    real(WP) :: m_snap

    if (PRESENT(m)) then
       $CHECK_BOUNDS(SIZE(m),SIZE(x))
    endif

    ! Snap model points to fix possible numerical issues

    ! Central point

    if (x(1) > 0._WP .AND. x(1) < dx_snap) then

       x(1) = 0._WP

       if (PRESENT(m)) then
          m(1) = 0._WP
       endif

       if (check_log_level('INFO')) then
          write(OUTPUT_UNIT, 100) 'Snapping central point to x=0'
100       format(3X,A)
       endif

    endif

    ! Other points

    snap_loop : do i = 2, SIZE(x)-1

       if (x(i+1) - x(i) > 0._WP .AND. x(i+1) - x(i) < dx_snap) then
          
          x_snap = 0.5_WP*(x(i+1) + x(i))
          x(i:i+1) = x_snap

          if (PRESENT(m)) then
             m_snap = 0.5_WP*(m(i+1) + m(i))
             m(i:i+1) = m_snap
          endif

          if (check_log_level('INFO')) then
             write(OUTPUT_UNIT, 110) 'Snapping points', i, 'and', i+1, 'to x=', x_snap
110          format(3X,A,1X,I0,1X,A,1X,I0,1X,A,F6.4)
          endif
             
       end if

    end do snap_loop

    ! Finish

    return

  end subroutine snap_points

end module gyre_model_util
