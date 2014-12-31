! Program  : gyre_input
! Purpose  : input routines
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

module gyre_input

  ! Uses

  use core_kinds
  use core_order
  use core_parallel
  use core_system

  use gyre_constants

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: init_system
  public :: read_constants
  public :: read_model

contains

  subroutine init_system (filename, gyre_dir)

    character(:), allocatable, intent(out) :: filename
    character(:), allocatable, intent(out) :: gyre_dir

    integer :: status

    ! Get command-line arguments

    $ASSERT(n_arg() == 1,Invalid number of arguments)

    call get_arg(1, filename)

    ! Get environment variables

    call get_env('GYRE_DIR', gyre_dir, status)
    $ASSERT(status == 0,The GYRE_DIR environment variable is not set)

    ! Finish

    return

  end subroutine init_system

!****

  subroutine read_model (unit, x_bc, ml)

    use gyre_model
    use gyre_evol_model
    use gyre_scons_model
    use gyre_poly_model
    use gyre_hom_model
    use gyre_mesa_file
    use gyre_osc_file
    use gyre_losc_file
    use gyre_fgong_file
    use gyre_famdl_file
    use gyre_amdl_file
    $if ($HDF5)
    use gyre_b3_file
    use gyre_gsm_file
    use gyre_poly_file
    $endif

    integer, intent(in)                  :: unit
    real(WP), allocatable, intent(out)   :: x_bc(:)
    class(model_t), pointer, intent(out) :: ml

    integer                 :: n_ml
    character(256)          :: model_type
    character(256)          :: file_format
    character(256)          :: data_format
    character(256)          :: deriv_type
    character(FILENAME_LEN) :: file
    real(WP)                :: Gamma_1
    logical                 :: reconstruct_As
    logical                 :: uniform_rot
    real(WP)                :: Omega_rot
    type(evol_model_t)      :: ec
    type(scons_model_t)     :: sc
    type(poly_model_t)      :: pc
    type(hom_model_t)       :: hc

    namelist /model/ model_type, file_format, data_format, deriv_type, &
                     file, Gamma_1, &
                     reconstruct_As, uniform_rot, Omega_rot

    ! Count the number of model namelists

    rewind(unit)

    n_ml = 0

    count_loop : do
       read(unit, NML=model, END=100)
       n_ml = n_ml + 1
    end do count_loop

100 continue

    $ASSERT(n_ml == 1,Input file should contain exactly one &model namelist)

    ! Read model parameters

    model_type = ''
    file_format = ''
    data_format = ''
    deriv_type = 'MONO'
    uniform_rot = .FALSE.
    reconstruct_As = .FALSE.

    file = ''

    Gamma_1 = 5._WP/3._WP
    Omega_rot = 0._WP

    rewind(unit)
    read(unit, NML=model)

    ! Read/initialize the model

    select case (model_type)
    case ('EVOL')

       select case (file_format)
       case ('MESA')
          if (uniform_rot) then
             call read_mesa_model(file, deriv_type, ec, x=x_bc, uni_Omega_rot=Omega_rot)
          else
             call read_mesa_model(file, deriv_type, ec, x=x_bc)
          endif
       case('B3')
          $if($HDF5)
          call read_b3_model(file, deriv_type, ec, x=x_bc)
          $else
          $ABORT(No HDF5 support, therefore cannot read B3-format files)
          $endif
       case ('GSM')
          $if($HDF5)
          call read_gsm_model(file, deriv_type, ec, x=x_bc)
          $else
          $ABORT(No HDF5 support, therefore cannot read GSM-format files)
          $endif
       case ('OSC')
          call read_osc_model(file, deriv_type, data_format, ec, x=x_bc)
       case ('LOSC')
          call read_losc_model(file, deriv_type, ec, x=x_bc)
       case ('FGONG')
          call read_fgong_model(file, deriv_type, data_format, ec, x=x_bc) 
       case ('FAMDL')
          call read_famdl_model(file, deriv_type, data_format, ec, x=x_bc)
       case ('AMDL')
          call read_amdl_model(file, deriv_type, ec, x=x_bc)
       case default
          $ABORT(Invalid file_format)
       end select

       ec%reconstruct_As = reconstruct_As

       allocate(ml, SOURCE=ec)

    case ('SCONS')

       select case (file_format)
       case ('MESA')
          call read_mesa_model(file, sc, x=x_bc)
       case ('FGONG')
          call read_fgong_model(file, data_format, sc, x=x_bc)
       case default
          $ABORT(Invalid file_format)
       end select

       allocate(ml, SOURCE=sc)

    case ('POLY')

       $if($HDF5)
       call read_poly_model(file, deriv_type, pc, x=x_bc)
       $else
       $ABORT(No HDF5 support, therefore cannot read POLY files)
       $endif

       allocate(ml, SOURCE=pc)

    case ('HOM')

       hc = hom_model_t(Gamma_1, Omega_rot)

       allocate(ml, SOURCE=hc)

    case default

       $ABORT(Invalid model_type)

    end select

    ! Finish

    return

  end subroutine read_model

!****

  subroutine read_constants (unit)

    integer, intent(in) :: unit

    integer :: n_cn

    namelist /constants/ G_GRAVITY, C_LIGHT, A_RADIATION, &
                         M_SUN, R_SUN, L_SUN

    ! Count the number of constants namelists

    rewind(unit)

    n_cn = 0

    count_loop : do
       read(unit, NML=constants, END=100)
       n_cn = n_cn + 1
    end do count_loop

100 continue

    $ASSERT(n_cn == 1,Input file should contain exactly one &constants namelist)

    ! Read constants

    rewind(unit)
    read(unit, NML=constants)

    ! Finish

    return

  end subroutine read_constants

end module gyre_input
