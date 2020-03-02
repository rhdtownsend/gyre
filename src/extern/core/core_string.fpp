! Module   : core_string
! Purpose  : string manipulation

$include 'core.inc'
$include 'core_parallel.inc'

module core_string

  ! Uses

  use core_kinds

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameters

  integer, parameter :: GET_BUFFER_LEN = 256

  ! Interfaces

  interface get
     module procedure get_
     module procedure get_unit_
     module procedure get_set_
     module procedure get_unit_set_
  end interface get

  interface put
     module procedure put_
     module procedure put_unit_
  end interface put
     
  interface put_line
     module procedure put_line_
     module procedure put_line_unit_
  end interface put_line

  interface replace
     module procedure replace_auto_
     module procedure replace_fixed_
     module procedure replace_target_
  end interface replace

  ! Access specifiers

  private

  public :: get
  public :: put
  public :: extract
  public :: insert
  public :: remove
  public :: replace
  public :: split

  ! Procedures

contains

  subroutine get_ (string, maxlen, iostat)

    character(:), allocatable, intent(out) :: string
    integer, optional, intent(in)          :: maxlen
    integer, optional, intent(out)         :: iostat

    ! Read from the default unit into a string

    call get(INPUT_UNIT, string, maxlen, iostat)

    ! Finish

    return

  end subroutine get_

!****

  subroutine get_unit_ (unit, string, maxlen, iostat)

    integer, intent(in)                    :: unit
    character(:), allocatable, intent(out) :: string
    integer, optional, intent(in)          :: maxlen
    integer, optional, intent(out)         :: iostat

    integer                   :: n_chars_remain
    integer                   :: n_chars_read
    character(GET_BUFFER_LEN) :: buffer
    integer                   :: local_iostat

    if (PRESENT(maxlen)) then
       n_chars_remain = maxlen
    else
       n_chars_remain = HUGE(0)
    endif

    ! Read from the specified unit into a string

    string = ''

    read_loop : do

       if (n_chars_remain <= 0) return

       n_chars_read = MIN(n_chars_remain, GET_BUFFER_LEN)

       if (PRESENT(iostat)) then
          read(unit=unit, FMT='(A)', ADVANCE='NO', &
               IOSTAT=iostat, SIZE=n_chars_read) buffer(:n_chars_read)
          if (IS_IOSTAT_EOR(iostat)) then
             iostat = 0
             exit read_loop
          elseif (IS_IOSTAT_END(iostat)) then
             exit read_loop
          elseif (iostat > 0) then
             return
          endif
       else
          read(unit=unit, FMT='(A)', ADVANCE='NO', &
               IOSTAT=local_iostat, SIZE=n_chars_read) buffer(:n_chars_read)
          if (IS_IOSTAT_EOR(iostat) .OR. IS_IOSTAT_END(iostat)) exit read_loop
       endif

       string = string//buffer(:n_chars_read)
       n_chars_remain = n_chars_remain - n_chars_read

    end do read_loop

    string = string//buffer(:n_chars_read)

    ! Finish

    return

  end subroutine get_unit_

!****

  subroutine get_set_ (string, set, separator, maxlen, iostat)

    character(:), allocatable, intent(out)           :: string
    character(*), intent(in)                         :: set
    character(:), allocatable, optional, intent(out) :: separator
    integer, optional, intent(in)                    :: maxlen
    integer, optional, intent(out)                   :: iostat

    ! Read from the default unit into a string, with a custom
    ! character-string separator

    call get(INPUT_UNIT, string, set, separator, maxlen, iostat)

    ! Finish

    return

  end subroutine get_set_

!****

  subroutine get_unit_set_ (unit, string, set, separator, maxlen, iostat)

    integer, intent(in)                              :: unit
    character(:), allocatable, intent(out)           :: string
    character(*), intent(in)                         :: set
    character(:), allocatable, optional, intent(out) :: separator
    integer, optional, intent(in)                    :: maxlen
    integer, optional, intent(out)                   :: iostat

    integer      :: n_chars_remain
    character(1) :: buffer
    integer      :: i_set
    integer      :: local_iostat

    if (PRESENT(maxlen)) then
       n_chars_remain = maxlen
    else
       n_chars_remain = HUGE(1)
    endif

    if (PRESENT(separator)) separator = ''

    ! Read from the specified unit into a string, with a custom
    ! character-string separator

    string = ''

    read_loop : do

       if (n_chars_remain <= 0) return

       if (PRESENT(iostat)) then
          read(unit=unit, FMT='(A1)', ADVANCE='NO', IOSTAT=iostat) buffer
          if (iostat /= 0) exit read_loop
       else
          read(unit=unit, FMT='(A1)', ADVANCE='NO', IOSTAT=local_iostat) buffer
          if (local_iostat /= 0) exit read_loop
       endif

       i_set = SCAN(buffer, set)

       if (i_set == 1) then
          if (PRESENT(separator)) separator = buffer
          exit read_loop
       endif

       string = string//buffer
       n_chars_remain = n_chars_remain - 1

    end do read_loop

    ! Finish

    return

  end subroutine get_unit_set_

!****

  subroutine put_ (string, iostat)

    character(*), intent(in)       :: string
    integer, optional, intent(out) :: iostat

    ! Append a string to the current record of the default unit

    call put(OUTPUT_UNIT, string, iostat)

    ! Finish

  end subroutine put_

!****

  subroutine put_unit_ (unit, string, iostat)

    integer, intent(in)            :: unit
    character(*), intent(in)       :: string
    integer, optional, intent(out) :: iostat

    ! Append a string to the current record of the specified unit

    if (PRESENT(iostat)) then
       write(unit=unit, FMT='(A)', ADVANCE='NO', IOSTAT=iostat) string
    else
       write(unit=unit, FMT='(A)', ADVANCE='NO') string
    endif

    ! Finish

    return

  end subroutine put_unit_

!****

  subroutine put_line_ (string, iostat)

    character(*), intent(in)       :: string
    integer, optional, intent(out) :: iostat

    ! Append a string to the current record of the default unit,
    ! terminating the record

    call put_line(OUTPUT_UNIT, string, iostat)

    ! Finish

    return

  end subroutine put_line_

!****

  subroutine put_line_unit_ (unit, string, iostat)

    integer, intent(in)            :: unit
    character(*), intent(in)       :: string
    integer, optional, intent(out) :: iostat

    ! Append a string to the current record of the specified unit,
    ! terminating the record

    if (PRESENT(iostat)) then
       write(unit=unit, FMT='(A,/)', ADVANCE='NO', IOSTAT=iostat) string
    else
       write(unit=unit, FMT='(A,/)', ADVANCE='NO') string
    endif

    ! Finish

    return

  end subroutine put_line_unit_

!****

  function extract (string, start, finish) result (substring)

    character(*), intent(in)      :: string
    integer, optional, intent(in) :: start
    integer, optional, intent(in) :: finish
    character(:), allocatable     :: substring

    integer :: start_
    integer :: finish_

    if (PRESENT(start)) then
       start_ = MAX(1, start)
    else
       start_ = 1
    endif

    if (PRESENT(finish)) then
       finish_ = MIN(LEN(string), finish)
    else
       finish_ = LEN(string)
    endif

    ! Extract a substring from a string

    substring = string(start_:finish_)

    ! Finish

    return

  end function extract

!****

  function insert (string, start, substring) result (new_string)

    character(*), intent(in)  :: string
    integer, intent(in)       :: start
    character(*), intent(in)  :: substring
    character(:), allocatable :: new_string

    integer :: start_

    ! Insert a substring into a string

    start_ = MAX(1, MIN(start, LEN(string)+1))

    new_string = string(:start_-1)//substring//string(start_:)

    ! Finish

    return

  end function insert

!****

  function remove (string, start, finish) result (new_string)

    character(*), intent(in)      :: string
    integer, optional, intent(in) :: start
    integer, optional, intent(in) :: finish
    character(:), allocatable     :: new_string

    integer :: start_
    integer :: finish_

    if (PRESENT(start)) then
       start_ = MAX(1, start)
    else
       start_ = 1
    endif

    if (PRESENT(finish)) then
       finish_ = MIN(LEN(string), finish)
    else
       finish_ = LEN(string)
    endif

    ! Remove a substring from a string

    if (finish_ >= start_) then
       new_string = string(:start_-1)//string(finish_+1:)
    else
       new_string = string
    endif

    ! Finish

    return

  end function remove

!****

  function replace_auto_ (string, start, substring) result (new_string)

    character(*), intent(in)  :: string
    integer, intent(in)       :: start
    character(*), intent(in)  :: substring
    character(:), allocatable :: new_string

    ! Replace part of a string with a substring

    new_string = replace(string, start, MAX(start, 1)+LEN(substring)-1, substring)

    ! Finish

    return

  end function replace_auto_

!****

  function replace_fixed_ (string, start, finish, substring) result (new_string)

    character(*), intent(in)  :: string
    integer, intent(in)       :: start
    integer, intent(in)       :: finish
    character(*), intent(in)  :: substring
    character(:), allocatable :: new_string

    integer :: start_
    integer :: finish_

    ! Replace part of a string with a substring

    start_ = MAX(1, start)
    finish_ = MIN(LEN(string), finish)

    if (finish_ < start_) then
       new_string = insert(string, start_, substring)
    else
       new_string = string(:start_-1)//substring//string(finish_+1:)
    endif

    ! Finish

    return

  end function replace_fixed_

!****

  function replace_target_ (string, target, substring, every, back) result (new_string)

    character(*), intent(in)      :: string
    character(*), intent(in)      :: target
    character(*), intent(in)      :: substring
    logical, optional, intent(in) :: every
    logical, optional, intent(in) :: back
    character(:), allocatable     :: new_string

    logical                   :: every_
    logical                   :: back_
    character(:), allocatable :: work_string
    integer                   :: length_target
    integer                   :: i_target

    if (PRESENT(every)) then
       every_ = every
    else
       every_ = .FALSE.
    endif

    if (PRESENT(back)) then
       back_ = back
    else
       back_ = .FALSE.
    endif

    $ASSERT(LEN(target) /= 0,Invalid target length)

    ! Replace part of a string with a substring, at a location
    ! matching a string target

    new_string = ''

    work_string = string

    length_target = LEN(target)

    replace_loop : do

       i_target = INDEX(work_string, target, back_)

       if (i_target == 0) exit replace_loop

       if (back_) then
          new_string = substring//extract(work_string, start=i_target+length_target)//new_string
          work_string = extract(work_string, finish=i_target-1)
       else
          new_string = new_string//extract(work_string, finish=i_target-1)//substring
          work_string = extract(work_string, start=i_target+length_target)
       endif

       if (.NOT. every_) exit replace_loop

    end do replace_loop

    if (back_) then
       new_string = work_string//new_string
    else
       new_string = new_string//work_string
    endif

    ! Finish

    return

  end function replace_target_

!****

  subroutine split (string, word, set, separator, back)

    character(:), allocatable, intent(inout)         :: string
    character(:), allocatable, intent(out)           :: word
    character(*), intent(in)                         :: set
    character(:), allocatable, optional, intent(out) :: separator
    logical, optional, intent(in)                    :: back

    logical :: back_
    integer :: i_separator

    if (PRESENT(back)) then
       back_ = back
    else
       back_ = .FALSE.
    endif

    ! Split a string into two strings

    i_separator = SCAN(string, set, back_)

    if (i_separator /= 0) then

       if (back_) then
          word = extract(string, start=i_separator+1)
          if (PRESENT(separator)) separator = extract(string, start=i_separator, finish=i_separator)
          string = extract(string, finish=i_separator-1)
       else
          word = extract(string, finish=i_separator-1)
          if (PRESENT(separator)) separator = extract(string, start=i_separator, finish=i_separator)
          string = extract(string, start=i_separator+1)
       endif

    else

       word = string
       if (PRESENT(separator)) separator = ''
       string = ''

    endif

    ! Finish

    return

  end subroutine split

end module core_string
