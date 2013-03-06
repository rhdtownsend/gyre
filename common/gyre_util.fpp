! Module   : gyre_util
! Purpose  : misc utility routines

$include 'core.inc'

module gyre_util

  ! Uses

  use core_kinds
  use core_parallel

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: open_input
  public :: write_header

  ! Procedures

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

end module gyre_util
