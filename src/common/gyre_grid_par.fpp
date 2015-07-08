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
     real(WP)                :: x_i = 0._WP
     real(WP)                :: x_o = 1._WP
     real(WP)                :: alpha_osc = 0._WP
     real(WP)                :: alpha_exp = 0._WP
     real(WP)                :: alpha_thm = 0._WP
     real(WP)                :: alpha_str = 0._WP
     real(WP)                :: s = 0
     integer                 :: n = 0
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
  public :: read_shoot_grid_par
  public :: read_recon_grid_par
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

!****

  $define $READ_GRID_PAR $sub

  $local $NAME $1

  subroutine read_${NAME}_grid_par (unit, gp)

    integer, intent(in)                        :: unit
    type(grid_par_t), allocatable, intent(out) :: gp(:)

    integer                     :: n_gp
    integer                     :: i
    real(WP)                    :: x_i
    real(WP)                    :: x_o
    real(WP)                    :: alpha_osc
    real(WP)                    :: alpha_exp
    real(WP)                    :: alpha_thm
    real(WP)                    :: alpha_str
    real(WP)                    :: s
    integer                     :: n
    character(LEN(gp%file))     :: file
    character(LEN(gp%op_type))  :: op_type
    character(LEN(gp%tag_list)) :: tag_list

    namelist /${NAME}_grid/ x_i, x_o, &
                            alpha_osc, alpha_exp, alpha_thm, alpha_str, &
                            s, n, file, op_type, tag_list

    ! Count the number of grid namelists

    rewind(unit)

    n_gp = 0

    count_loop : do
       read(unit, NML=${NAME}_grid, END=100)
       n_gp = n_gp + 1
    end do count_loop

100 continue

    ! Read grid parameters

    rewind(unit)

    allocate(gp(n_gp))

    read_loop : do i = 1, n_gp

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

       gp(i) = grid_par_t(x_i=x_i, &
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

  end subroutine read_${NAME}_grid_par

  $endsub

  $READ_GRID_PAR(shoot)
  $READ_GRID_PAR(recon)

!****

  $if ($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_${RANK}_ (gp, root_rank)

    type(grid_par_t), intent(inout) :: gp$ARRAY_SPEC($RANK)
    integer, intent(in)             :: root_rank

    ! Broadcast the grid_par_t

    call bcast(gp%x_i, root_rank)
    call bcast(gp%x_o, root_rank)

    call bcast(gp%alpha_osc, root_rank)
    call bcast(gp%alpha_exp, root_rank)
    call bcast(gp%alpha_thm, root_rank)
    call bcast(gp%alpha_str, root_rank)

    call bcast(gp%s, root_rank)

    call bcast(gp%n, root_rank)

    call bcast(gp%file, root_rank)

    call bcast(gp%op_type, root_rank)
    call bcast(gp%tag_list, root_rank)

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
