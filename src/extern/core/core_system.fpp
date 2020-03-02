! Module   : core_system
! Purpose  : operating system support

$include 'core.inc'

module core_system

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Interfaces

  interface get_arg
     module procedure get_arg_i_i4_
     module procedure get_arg_i_i8_
     module procedure get_arg_r_sp_
     module procedure get_arg_r_dp_
     module procedure get_arg_c_sp_
     module procedure get_arg_c_dp_
     module procedure get_arg_a_
     module procedure get_arg_l_
  end interface get_arg

  ! Access specifiers

  private

  public :: n_arg
  public :: get_arg

contains

  function n_arg ()

    integer :: n_arg

    ! Get the number of arguments

    n_arg = COMMAND_ARGUMENT_COUNT()

    ! Finish

    return

  end function n_arg

!****

  $define $GET_ARG $sub

  $local $INFIX $1
  $local $ARG_TYPE $2

  subroutine get_arg_${INFIX}_ (i_arg, arg)

    integer, intent(in)    :: i_arg
    $ARG_TYPE, intent(out) :: arg

    integer                       :: buffer_len
    character(LEN=:), allocatable :: buffer

    ! Read the argument into a character buffer

    call GET_COMMAND_ARGUMENT(i_arg, LENGTH=buffer_len)

    allocate(character(LEN=buffer_len)::buffer)

    call GET_COMMAND_ARGUMENT(i_arg, buffer)

    ! Extract the argument from the buffer

    read(buffer, *) arg

    ! Finish

    return

  end subroutine get_arg_${INFIX}_

  $endsub

  $GET_ARG(i_i4,integer(I4))
  $GET_ARG(i_i8,integer(I8))
  $GET_ARG(r_sp,real(SP))
  $GET_ARG(r_dp,real(DP))
  $GET_ARG(c_sp,complex(SP))
  $GET_ARG(c_dp,complex(DP))
  $GET_ARG(a,character(LEN=*))
  $GET_ARG(l,logical)
  
end module core_system
