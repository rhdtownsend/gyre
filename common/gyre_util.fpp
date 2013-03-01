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

  public :: write_header

  ! Procedures

contains

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
