! Module   : core_env
! Purpose  : environment access

$include 'core.inc'

module core_env

  ! Uses

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: get_env

  ! Procedures

contains

  subroutine get_env (name, value, trim_name)

    character(*), intent(in)               :: name
    character(:), allocatable, intent(out) :: value
    logical, optional, intent(in)          :: trim_name

    integer :: length

    ! Get the named environment variable

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, TRIM_NAME=trim_name)

    if (length > 0) then

       allocate(character(length)::value)

       call GET_ENVIRONMENT_VARIABLE(name, VALUE=value, TRIM_NAME=trim_name)

    endif

    ! Finish

    return

  end subroutine get_env

end module core_env
