! Module   : gyre_grid_par
! Purpose  : grid parameters
!
! Copyright 2013-2020 Rich Townsend & The GYRE Team
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
     real(WP)        :: alpha_ctr = 0._WP
     real(WP)        :: alpha_thm = 0._WP
     real(WP)        :: alpha_str = 0._WP
     real(WP)        :: dx_min = sqrt(EPSILON(0._WP))
     real(WP)        :: dx_max = HUGE(0._WP)
     integer         :: n_iter_max = 32
     logical         :: resolve_ctr = .TRUE.
     character (64)  :: isolation = 'NONE'
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
    real(WP)                      :: alpha_ctr
    real(WP)                      :: alpha_thm
    real(WP)                      :: alpha_str
    real(WP)                      :: dx_min
    real(WP)                      :: dx_max
    integer                       :: n_iter_max
    logical                       :: resolve_ctr
    character(LEN(gr_p%isolation)):: isolation
    character(LEN(gr_p%tag_list)) :: tag_list

    namelist /grid/ x_i, x_o, alpha_osc, alpha_exp, alpha_ctr, &
         alpha_thm, alpha_str, dx_min, dx_max, n_iter_max, &
         resolve_ctr, isolation, tag_list

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
       alpha_ctr = gr_p(i)%alpha_ctr
       alpha_thm = gr_p(i)%alpha_thm
       alpha_str = gr_p(i)%alpha_str
       dx_min = gr_p(i)%dx_min
       dx_max = gr_p(i)%dx_max
       n_iter_max = gr_p(i)%n_iter_max
       isolation = gr_p(i)%isolation
       resolve_ctr = gr_p(i)%resolve_ctr
       tag_list = gr_p(i)%tag_list

       ! Read the namelist

       read(unit, NML=grid)

       ! Store read values

       gr_p(i)%x_i = x_i
       gr_p(i)%x_o = x_o
       gr_p(i)%alpha_osc = alpha_osc
       gr_p(i)%alpha_exp = alpha_exp
       gr_p(i)%alpha_ctr = alpha_ctr
       gr_p(i)%alpha_thm = alpha_thm
       gr_p(i)%alpha_str = alpha_str
       gr_p(i)%dx_min = dx_min
       gr_p(i)%dx_max = dx_max
       gr_p(i)%n_iter_max = n_iter_max
       gr_p(i)%isolation = isolation
       gr_p(i)%resolve_ctr = resolve_ctr
       gr_p(i)%tag_list = tag_list

    end do read_loop

    ! Finish

    return

  end subroutine read_grid_par

end module gyre_grid_par
