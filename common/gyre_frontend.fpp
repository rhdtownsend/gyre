! Program  : gyre_frontend
! Purpose  : frontend routines
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

module gyre_frontend

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel

  use gyre_oscpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: open_input
  public :: write_header
  public :: init_coeffs
  public :: init_oscpar

contains

  subroutine open_input (unit)

    integer, intent(out) :: unit

    character(LEN=1024) :: line

    ! Make standard input available through a scratch file

    open(NEWUNIT=unit, STATUS='SCRATCH')

    do
       read(INPUT_UNIT, 100, END=200) line
100    format(A)
       write(unit, *) line
    end do

200 continue

    ! Finish

    return

  end subroutine open_input

!****

  subroutine write_header (header, underchar)

    character(LEN=*), intent(in)           :: header
    character(LEN=*), intent(in), optional :: underchar

    ! Write out the header

    if(MPI_RANK == 0) then

       write(OUTPUT_UNIT, '()')

       write(OUTPUT_UNIT, '(A)') header

       if(PRESENT(underchar)) then
          if(underchar == '') then
             write(OUTPUT_UNIT, '(A)') REPEAT(' ', LEN_TRIM(header))
          else
             write(OUTPUT_UNIT, '(A)') REPEAT(underchar, LEN_TRIM(header)/LEN_TRIM(underchar))
          endif
       endif

       write(OUTPUT_UNIT, '()')

    endif

    ! Finish

    return

  end subroutine write_header

!****

  subroutine init_coeffs (unit, x_mc, mc, tc)

    use gyre_mech_coeffs
    use gyre_therm_coeffs
    use gyre_mech_coeffs_mpi
    use gyre_therm_coeffs_mpi
    use gyre_mesa_file
    use gyre_fgong_file
    use gyre_osc_file
    use gyre_b3_file
    use gyre_mech_coeffs_poly
    use gyre_mech_coeffs_hom

    integer, intent(in)                                         :: unit
    real(WP), allocatable, intent(out)                          :: x_mc(:)
    $if($GFORTRAN_PR56218)
    class(mech_coeffs_t), allocatable, intent(inout)            :: mc
    class(therm_coeffs_t), allocatable, intent(inout), optional :: tc
    $else
    class(mech_coeffs_t), allocatable, intent(out)              :: mc
    class(therm_coeffs_t), allocatable, intent(out), optional   :: tc
    $endif

    character(LEN=256)          :: coeffs_type
    character(LEN=FILENAME_LEN) :: file
    real(WP)                    :: G
    real(WP)                    :: Gamma_1

    namelist /coeffs/ coeffs_type, file, G, Gamma_1

    ! Read structure coefficients parameters

    if(MPI_RANK == 0) then

       coeffs_type = ''

       file = ''

       G = G_GRAVITY
       Gamma_1 = 5._WP/3._WP

       rewind(unit)
       read(unit, NML=coeffs)

    endif

    $if($MPI)
    call bcast(coeffs_type, 0)
    $endif

    ! Read/initialize the mech_coeffs

    if(MPI_RANK == 0) then

       select case(coeffs_type)
       case('MESA')
          call read_mesa_file(file, G, mc, tc, x=x_mc)
       case('B3')
          call read_b3_file(file, G, mc, tc, x=x_mc)
       case('FGONG')
          call read_fgong_file(file, G, mc, x=x_mc) 
       case('OSC')
          call read_osc_file(file, G, mc, x=x_mc)
       case('POLY')
          allocate(mech_coeffs_poly_t::mc)
          select type (mc)
          type is (mech_coeffs_poly_t)
             call mc%read(file, x_mc)
          end select
       case('HOM')
          allocate(mech_coeffs_hom_t::mc)
          select type (mc)
          type is (mech_coeffs_hom_t)
             call mc%init(Gamma_1)
          end select
       end select

    endif

    $if($MPI)
    call alloc_bcast(mc, 0)
    if(PRESENT(tc)) call alloc_bcast(tc, 0)
    call alloc_bcast(x_mc, 0)
    $endif

    ! Finish

    return

  end subroutine init_coeffs

!****

  subroutine init_oscpar (unit, op)

    integer, intent(in)         :: unit
    type(oscpar_t), intent(out) :: op

    integer            :: l
    character(LEN=256) :: outer_bound_type

    namelist /oscpar/ l, outer_bound_type

    ! Read oscillation parameters

    if(MPI_RANK == 0) then

       l = 0

       outer_bound_type = 'ZERO'

       rewind(unit)
       read(unit, NML=oscpar)

    endif

    $if($MPI)
    call bcast(l, 0)
    call bcast(outer_bound_type, 0)
    $endif

    ! Initialize the oscilaltion parameters

    call op%init(l, outer_bound_type)

    ! Finish

    return

  end subroutine init_oscpar

end module gyre_frontend
