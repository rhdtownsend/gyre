! Module   : gyre_grid_par
! Purpose  : grid parameters
!
! Copyright 2013-2015 Rich Townsend
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
$include 'core_parallel.inc'

module gyre_grid_par

  ! Uses

  use core_kinds
  use gyre_constants
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: grid_par_t
     real(WP)                :: x_i
     real(WP)                :: x_o
     real(WP)                :: alpha_osc
     real(WP)                :: alpha_exp
     real(WP)                :: alpha_thm
     real(WP)                :: alpha_str
     real(WP)                :: s
     integer                 :: n
     character(FILENAME_LEN) :: file
     character(64)           :: op_type
     character(2048)         :: tag_list
  end type grid_par_t

  ! Interfaces

  $if ($MPI)

  interface bcast
     module procedure bcast_0_
     module procedure bcast_1_
  end interface bcast

  interface bcast_alloc
     module procedure bcast_alloc_0_
     module procedure bcast_alloc_1_
  end interface bcast_alloc

  $endif

  ! Access specifiers

  private

  public :: grid_par_t
  public :: read_grid_par
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  subroutine read_grid_par (unit, gd_p)

    integer, intent(in)                        :: unit
    type(grid_par_t), allocatable, intent(out) :: gd_p(:)

    integer                       :: n_gd_p
    integer                       :: i
    real(WP)                      :: x_i
    real(WP)                      :: x_o
    real(WP)                      :: alpha_osc
    real(WP)                      :: alpha_exp
    real(WP)                      :: alpha_thm
    real(WP)                      :: alpha_str
    real(WP)                      :: s
    integer                       :: n
    character(LEN(gd_p%file))     :: file
    character(LEN(gd_p%op_type))  :: op_type
    character(LEN(gd_p%tag_list)) :: tag_list

    namelist /${NAME}_grid/ x_i, x_o, &
                            alpha_osc, alpha_exp, alpha_thm, alpha_str, &
                            s, n, file, op_type, tag_list

    ! Count the number of grid namelists

    rewind(unit)

    n_gd_p = 0

    count_loop : do
       read(unit, NML=${NAME}_grid, END=100)
       n_gd_p = n_gd_p + 1
    end do count_loop

100 continue

    ! Read grid parameters

    rewind(unit)

    allocate(gd_p(n_gd_p))

    read_loop : do i = 1, n_gd_p

       x_i = 0._WP
       x_o = 1._WP

       alpha_osc = 0._WP
       alpha_exp = 0._WP
       alpha_thm = 0._WP
       alpha_str = 0._WP

       s = 0._WP

       n = 0

       file = ''

       op_type = 'CREATE_CLONE'
       tag_list = ''

       read(unit, NML=${NAME}_grid)

       ! Initialize the grid_par

       gd_p(i) = grid_par_t(x_i=x_i, &
                            x_o=x_o, &
                            alpha_osc=alpha_osc, &
                            alpha_exp=alpha_exp, &
                            alpha_thm=alpha_thm, alpha_str=alpha_str, &
                            s=s, &
                            n=n, &
                            file=file, &
                            op_type=op_type, &
                            tag_list=tag_list)

    end do read_loop

    ! Finish

    return

  end subroutine read_grid_par

!****

  $if ($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_${RANK}_ (gd_p, root_rank)

    type(grid_par_t), intent(inout) :: gd_p$ARRAY_SPEC($RANK)
    integer, intent(in)             :: root_rank

    ! Broadcast the grid_par_t

    call bcast(gd_p%x_i, root_rank)
    call bcast(gd_p%x_o, root_rank)

    call bcast(gd_p%alpha_osc, root_rank)
    call bcast(gd_p%alpha_exp, root_rank)
    call bcast(gd_p%alpha_thm, root_rank)
    call bcast(gd_p%alpha_str, root_rank)

    call bcast(gd_p%s, root_rank)

    call bcast(gd_p%n, root_rank)

    call bcast(gd_p%file, root_rank)

    call bcast(gd_p%op_type, root_rank)
    call bcast(gd_p%tag_list, root_rank)

    ! Finish

    return

  end subroutine bcast_${RANK}_

  $endsub

  $BCAST(0)
  $BCAST(1)

!****

  $BCAST_ALLOC(type(grid_par_t),0)
  $BCAST_ALLOC(type(grid_par_t),1)

  $endif

end module gyre_grid_par
