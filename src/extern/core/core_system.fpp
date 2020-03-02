! Module   : core_system
! Purpose  : operating system support

$include 'core.inc'

module core_system

  ! Uses

  use core_kinds

  use ISO_FORTRAN_ENV

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

  interface get_env
     module procedure get_env_i_i4_
     module procedure get_env_i_i8_
     module procedure get_env_r_sp_
     module procedure get_env_r_dp_
     module procedure get_env_c_sp_
     module procedure get_env_c_dp_
     module procedure get_env_a_
     module procedure get_env_l_
  end interface get_env

  ! Access specifiers

  private

  public :: n_arg
  public :: get_arg
  public :: get_env

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
  $local $VALUE_TYPE $2

  subroutine get_arg_${INFIX}_ (number, value, status)

    integer, intent(in)            :: number
    $VALUE_TYPE, intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the numbered command argument

    ! Determine the argument length

    call GET_COMMAND_ARGUMENT(number, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else
       $ASSERT(status_ == 0,Error when reading command argument)
    endif

    ! Read the argument into a character buffer

    allocate(character(length)::buffer)

    call GET_COMMAND_ARGUMENT(number, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

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
  $GET_ARG(l,logical)

!****

  subroutine get_arg_a_ (number, value, status)

    integer, intent(in)                    :: number
    character(:), allocatable, intent(out) :: value
    integer, optional, intent(out)         :: status

    integer :: length
    integer :: status_

    ! Read the numbered command argument

    ! Determine the argument length

    call GET_COMMAND_ARGUMENT(number, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else
       $ASSERT(status_ == 0,Error when reading command argument)
    endif

    ! Read the argument

    if (length > 0) then
       allocate(character(length)::value)
       call GET_COMMAND_ARGUMENT(number, value)
    endif

    ! Finish

    return

  end subroutine get_arg_a_

!****

  $define $GET_ENV $sub

  $local $INFIX $1
  $local $VALUE_TYPE $2

  subroutine get_env_${INFIX}_ (name, value, status)

    character(*), intent(in)       :: name
    $VALUE_TYPE, intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the named environment variable

    ! Determine the variable length

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else
       $ASSERT(status_ == 0,Error when reading environment variable)
    endif

    ! Read the variable into a character buffer

    allocate(character(length)::buffer)

    call GET_ENVIRONMENT_VARIABLE(name, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_env_${INFIX}_

  $endsub

  $GET_ENV(i_i4,integer(I4))
  $GET_ENV(i_i8,integer(I8))
  $GET_ENV(r_sp,real(SP))
  $GET_ENV(r_dp,real(DP))
  $GET_ENV(c_sp,complex(SP))
  $GET_ENV(c_dp,complex(DP))
  $GET_ENV(l,logical)
  
!****

  subroutine get_env_a_ (name, value, status)

    character(*), intent(in)               :: name
    character(:), allocatable, intent(out) :: value
    integer, optional, intent(out)         :: status

    integer :: length
    integer :: status_

    ! Read the named environment variable

    ! Determine the variable length

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else
       $ASSERT(status_ == 0,Error when reading environment variable)
    endif

    ! Read the variable

    if (length > 0) then
       allocate(character(length)::value)
       call GET_ENVIRONMENT_VARIABLE(name, value)
    endif

    ! Finish

    return

  end subroutine get_env_a_

end module core_system
