! Module   : gyre_scan
! Purpose  : scan frequencies construction
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

module gyre_scan

  ! Uses

  use core_kinds
  use core_constants
  use core_func
  use core_order

  use gyre_base_coeffs
  use gyre_oscpar
  use gyre_gridpar
  use gyre_scanpar
  use gyre_grid
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: build_scan

contains

  subroutine build_scan (sp, bc, op, gp, x_in, omega)

    type(scanpar_t), intent(in)        :: sp(:)
    type(gridpar_t), intent(in)        :: gp(:)
    class(base_coeffs_t), intent(in)   :: bc
    type(oscpar_t), intent(in)         :: op
    real(WP), allocatable, intent(in)  :: x_in(:)
    real(WP), allocatable, intent(out) :: omega(:)

    real(WP) :: x_i
    real(WP) :: x_o
    integer  :: i
    real(WP) :: omega_min
    real(WP) :: omega_max
    integer  :: j

    ! Determine the grid range

    call grid_range(gp, bc, op, x_in, x_i, x_o)

    ! Loop through scanpars

    allocate(omega(0))

    sp_loop : do i = 1,SIZE(sp)

       ! Set up the frequency grid

       omega_min = sp(i)%freq_min/freq_scale(bc, op, x_o, sp(i)%freq_units)
       omega_max = sp(i)%freq_max/freq_scale(bc, op, x_o, sp(i)%freq_units)
       
       select case(sp(i)%grid_type)
       case('LINEAR')
          omega = [omega,(((sp(i)%n_freq-j)*omega_min + (j-1)*omega_max)/(sp(i)%n_freq-1), j=1,sp(i)%n_freq)]
       case('INVERSE')
          omega = [omega,((sp(i)%n_freq-1)/((sp(i)%n_freq-j)/omega_min + (j-1)/omega_max), j=1,sp(i)%n_freq)]
       case default
          $ABORT(Invalid grid_type)
       end select

    end do sp_loop

    ! Sort the frequencies

    omega = omega(sort_indices(omega))

    ! Finish

    return

  end subroutine build_scan

end module gyre_scan
