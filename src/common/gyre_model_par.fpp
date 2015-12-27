! Module   : gyre_osc_par
! Purpose  : model parameters
!
! Copyright 2015 Rich Townsend
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

module gyre_model_par

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: model_par_t
     real(WP)                :: Gamma_1
     real(WP)                :: Omega_rot
     character(256)          :: model_type
     character(256)          :: file_format
     character(256)          :: data_format
     character(256)          :: deriv_type
     character(FILENAME_LEN) :: file
     logical                 :: add_center
     logical                 :: uniform_rot
  end type model_par_t
   
  ! Interfaces

  $if ($MPI)

  interface bcast
     module procedure bcast_
  end interface bcast

  $endif

 ! Access specifiers

  private

  public :: model_par_t
  public :: read_model_par
  $if ($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine read_model_par (unit, ml_p)

    integer, intent(in)            :: unit
    type(model_par_t), intent(out) :: ml_p

    integer                         :: n_ml_p
    real(WP)                        :: Gamma_1
    real(WP)                        :: Omega_rot
    character(LEN(ml_p%model_type)) :: model_type
    character(LEN(ml_p%model_type)) :: file_format
    character(LEN(ml_p%model_type)) :: data_format
    character(LEN(ml_p%model_type)) :: deriv_type
    character(LEN(ml_p%file))       :: file
    logical                         :: add_center
    logical                         :: uniform_rot

    namelist /model/ Gamma_1, Omega_rot, model_type, file_format, data_format, &
                     deriv_type, file, add_center, uniform_rot
    
    ! Count the number of model namelists

    rewind(unit)

    n_ml_p = 0

    count_loop : do
       read(unit, NML=model, END=100)
       n_ml_p = n_ml_p + 1
    end do count_loop

100 continue

    $ASSERT(n_ml_p == 1,Input file should contain exactly one &model namelist)

    ! Read model parameters

    Gamma_1 = 5._WP/3._WP
    Omega_rot = 0._WP

    model_type = ''
    file_format = ''
    data_format = ''
    deriv_type = 'MONO'
    file = ''

    add_center = .TRUE.
    uniform_rot = .FALSE.

    rewind(unit)
    read(unit, NML=model)

    ! Initialize the model_par

    ml_p = model_par_t(Gamma_1=Gamma_1, &
                       Omega_rot=Omega_rot, &
                       model_type=model_type, &
                       file_format=file_format, &
                       data_format=data_format, &
                       deriv_type=deriv_type, &
                       file=file, &
                       add_center=add_center, &
                       uniform_rot=uniform_rot)

    ! Finish

    return

  end subroutine read_model_par

  !****

  $if ($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_${RANK}_ (ml_p, root_rank)

    type(model_par_t), intent(inout) :: ml_p$ARRAY_SPEC($RANK)
    integer, intent(in)              :: root_rank

    ! Broadcast the out_par_t

    call bcast(ml_p%Gamma_1, root_rank)
    call bcast(ml_p%Omega_rot, root_rank)
    call bcast(ml_p%model_type, root_rank)
    call bcast(ml_p%file_format, root_rank)
    call bcast(ml_p%data_format, root_rank)
    call bcast(ml_p%deriv_type, root_rank)
    call bcast(ml_p%file, root_rank)
    call bcast(ml_p%add_center, root_rank)
    call bcast(ml_p%uniform_rot, root_rank)

    ! Finish

    return

  end subroutine bcast_${RANK}_

  $endsub

  $BCAST(0)
  $BCAST(1)

  !****

  $BCAST_ALLOC(type(model_par_t),0)
  $BCAST_ALLOC(type(model_par_t),1)

  $endif

end module gyre_model_par
