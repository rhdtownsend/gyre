! Module   : gyre_grid_par
! Purpose  : grid parameters
!
! Copyright 2013-2018 Rich Townsend
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

module gyre_grid_par

  ! Uses

  use core_kinds

  use gyre_math

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: grid_par_t
     real(WP)        :: x_i = -HUGE(0._WP)
     real(WP)        :: x_o = HUGE(0._WP)
     real(WP)        :: alpha_osc = 0._WP
     real(WP)        :: alpha_exp = 0._WP
     real(WP)        :: alpha_thm = 0._WP
     real(WP)        :: alpha_str = 0._WP
     real(WP)        :: dx_min = sqrt(EPSILON(0._WP))
     integer         :: n_inner = 0
     integer         :: n_floor = 0
     integer         :: n_iter_max = 8
     character(2048) :: tag_list = ''
  end type grid_par_t

  ! Access specifiers

  private

  public :: grid_par_t
  public :: read_grid_par

  ! Procedures

contains

  subroutine read_grid_par (unit, gr_p)

    integer, intent(in)                        :: unit
    type(grid_par_t), allocatable, intent(out) :: gr_p(:)

    integer                       :: n_gr_p
    integer                       :: i
    real(WP)                      :: x_i
    real(WP)                      :: x_o
    real(WP)                      :: alpha_osc
    real(WP)                      :: alpha_exp
    real(WP)                      :: alpha_thm
    real(WP)                      :: alpha_str
    real(WP)                      :: dx_min
    integer                       :: n_inner
    integer                       :: n_floor
    integer                       :: n_iter_max
    character(LEN(gr_p%tag_list)) :: tag_list

    namelist /grid/ x_i, x_o, alpha_osc, alpha_exp, alpha_thm, alpha_str, &
                    dx_min, n_inner, n_floor, n_iter_max, tag_list

    ! Count the number of grid namelists

    rewind(unit)

    n_gr_p = 0

    count_loop : do
       read(unit, NML=grid, END=100)
       n_gr_p = n_gr_p + 1
    end do count_loop

100 continue

    ! Read grid parameters

    rewind(unit)

    allocate(gr_p(n_gr_p))

    read_loop : do i = 1, n_gr_p

       ! Set default values

       gr_p(i) = grid_par_t()

       x_i = gr_p(i)%x_i
       x_o = gr_p(i)%x_o
       alpha_osc = gr_p(i)%alpha_osc
       alpha_exp = gr_p(i)%alpha_exp
       alpha_thm = gr_p(i)%alpha_thm
       alpha_str = gr_p(i)%alpha_str
       dx_min = gr_p(i)%dx_min
       n_inner = gr_p(i)%n_inner
       n_floor = gr_p(i)%n_floor
       n_iter_max = gr_p(i)%n_iter_max
       tag_list = gr_p(i)%tag_list

       ! Read the namelist

       read(unit, NML=grid)

       ! Store read values

       gr_p(i)%x_i = x_i
       gr_p(i)%x_o = x_o
       gr_p(i)%alpha_osc = alpha_osc
       gr_p(i)%alpha_exp = alpha_exp
       gr_p(i)%alpha_thm = alpha_thm
       gr_p(i)%alpha_str = alpha_str
       gr_p(i)%dx_min = dx_min
       gr_p(i)%n_inner = n_inner
       gr_p(i)%n_floor = n_floor
       gr_p(i)%n_iter_max = n_iter_max
       gr_p(i)%tag_list = tag_list

    end do read_loop

    ! Finish

    return

  end subroutine read_grid_par

end module gyre_grid_par
