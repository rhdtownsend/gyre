! Module   : gyre_input
! Purpose  : input buffering

$include 'core.inc'

module gyre_input

  ! Uses

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: open_input

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

end module gyre_input
