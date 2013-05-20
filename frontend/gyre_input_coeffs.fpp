! Program  : gyre_input_coeffs
! Purpose  : coefficient input routines
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

module gyre_input_coeffs

  ! Uses

  use core_kinds
  use core_constants

  use gyre_base_coeffs
  use gyre_hom_base_coeffs
  use gyre_therm_coeffs
  use gyre_mesa_file
  use gyre_b3_file
  use gyre_gsm_file
  use gyre_fgong_file
  use gyre_osc_file
  use gyre_poly_file

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_coeffs

contains

  subroutine read_coeffs (unit, x_bc, bc, tc)

    integer, intent(in)                                         :: unit
    real(WP), allocatable, intent(out)                          :: x_bc(:)
    $if($GFORTRAN_PR56218)
    class(base_coeffs_t), allocatable, intent(inout)            :: bc
    class(therm_coeffs_t), allocatable, intent(inout), optional :: tc
    $else
    class(base_coeffs_t), allocatable, intent(out)              :: bc
    class(therm_coeffs_t), allocatable, intent(out), optional   :: tc
    $endif

    character(LEN=256)          :: coeffs_type
    character(LEN=256)          :: file_format
    character(LEN=256)          :: deriv_type
    character(LEN=FILENAME_LEN) :: file
    real(WP)                    :: G
    real(WP)                    :: Gamma_1

    namelist /coeffs/ coeffs_type, file_format, deriv_type, file, G, Gamma_1

    ! Read structure coefficients parameters

    coeffs_type = ''
    file_format = ''
    deriv_type = 'MONO'

    file = ''

    G = G_GRAVITY
    Gamma_1 = 5._WP/3._WP

    rewind(unit)
    read(unit, NML=coeffs, END=900)

    ! Read/initialize the base_coeffs

    select case (coeffs_type)
    case ('EVOL')
       select case (file_format)
       case ('MESA')
          call read_mesa_file(file, G, deriv_type, bc, tc, x=x_bc)
       case('B3')
          call read_b3_file(file, G, deriv_type, bc, tc, x=x_bc)
       case ('GSM')
          call read_gsm_file(file, G, deriv_type, bc, tc, x=x_bc)
       case ('OSC')
          call read_osc_file(file, G, deriv_type, bc, tc, x=x_bc)
       case ('FGONG')
          call read_fgong_file(file, G, deriv_type, bc, x=x_bc) 
       case default
          $ABORT(Invalid file_format)
       end select
    case ('POLY')
       call read_poly_file(file, deriv_type, bc, x=x_bc)
    case ('HOM')
       allocate(hom_base_coeffs_t::bc)
       select type (bc)
       type is (hom_base_coeffs_t)
          call bc%init(Gamma_1)
       end select
    case default
       $ABORT(Invalid coeffs_type)
    end select

    ! Finish

    return

    ! Jump-in point for end-of-file

900 continue

    $ABORT(No &coeffs namelist in input file)

  end subroutine read_coeffs

end module gyre_input_coeffs
