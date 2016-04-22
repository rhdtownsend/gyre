! Module   : gyre_model_par
! Purpose  : model parameters
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

module gyre_model_par

  ! Uses

  use core_kinds
  use core_constants, only : FILENAME_LEN

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: model_par_t
     real(WP)                :: Gamma_1
     real(WP)                :: Omega_rot
     real(WP)                :: dx_snap
     real(WP)                :: x_i
     real(WP)                :: x_o
     real(WP)                :: s
     character(256)          :: model_type
     character(256)          :: grid_type
     character(256)          :: file_format
     character(256)          :: data_format
     character(256)          :: deriv_type
     character(256)          :: Omega_units
     character(FILENAME_LEN) :: file
     integer                 :: n
     logical                 :: add_center
     logical                 :: repair_As
     logical                 :: uniform_rot
  end type model_par_t
   
 ! Access specifiers

  private

  public :: model_par_t
  public :: read_model_par

  ! Procedures

contains

  subroutine read_model_par (unit, ml_p)

    integer, intent(in)            :: unit
    type(model_par_t), intent(out) :: ml_p

    integer                          :: n_ml_p
    real(WP)                         :: Gamma_1
    real(WP)                         :: Omega_rot
    real(WP)                         :: dx_snap
    real(WP)                         :: x_i
    real(WP)                         :: x_o
    real(WP)                         :: s
    character(LEN(ml_p%model_type))  :: model_type
    character(LEN(ml_p%grid_type))   :: grid_type
    character(LEN(ml_p%file_format)) :: file_format
    character(LEN(ml_p%data_format)) :: data_format
    character(LEN(ml_p%deriv_type))  :: deriv_type
    character(LEN(ml_p%Omega_units)) :: Omega_units
    character(LEN(ml_p%file))        :: file
    integer                          :: n
    logical                          :: add_center
    logical                          :: repair_As
    logical                          :: uniform_rot

    namelist /model/ Gamma_1, Omega_rot, dx_snap, x_i, x_o, s, &
                     model_type, grid_type, file_format, data_format, deriv_type, &
                     Omega_units, file, n, add_center, repair_As, uniform_rot
    
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
    dx_snap = 0._WP

    x_i = 0._WP
    x_o = 1._WP
    s = 1._WP

    model_type = ''
    grid_type = 'UNIFORM'
    file_format = ''
    data_format = ''
    deriv_type = 'MONO'
    Omega_units = 'NONE'
    file = ''

    n = 10

    add_center = .TRUE.
    repair_As = .FALSE.
    uniform_rot = .FALSE.

    rewind(unit)
    read(unit, NML=model)

    ! Initialize the model_par

    ml_p = model_par_t(Gamma_1=Gamma_1, &
                       Omega_rot=Omega_rot, &
                       dx_snap=dx_snap, &
                       x_i=x_i, &
                       x_o=x_o, &
                       s=s, &
                       model_type=model_type, &
                       grid_type=grid_type, &
                       file_format=file_format, &
                       data_format=data_format, &
                       deriv_type=deriv_type, &
                       Omega_units=Omega_units, &
                       file=file, &
                       n=n, &
                       add_center=add_center, &
                       repair_As=repair_As, &
                       uniform_rot=uniform_rot)

    ! Finish

    return

  end subroutine read_model_par

end module gyre_model_par
