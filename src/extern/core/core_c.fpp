! Module   : core_c
! Purpose  : C interoperability

$include 'core.inc'

module core_c

  ! Uses

  use core_kinds

  use ISO_C_BINDING
  
  ! No implicit typing

  implicit none

  ! Interfaces

  interface
     function c_strlen (str) result (len) bind (C, name="strlen")
       use ISO_C_BINDING
       type(C_PTR), value :: str
       integer(C_SIZE_T)  :: len
     end function c_strlen
  end interface

  ! Access specifiers

  private

  public :: c_strlen
  public :: c_f_string

  ! Procedures

contains

  function c_f_string (c_str) result (f_str)
    
    type(C_PTR), intent(in)   :: c_str
    character(:), allocatable :: f_str

    integer                      :: l
    character(1,C_CHAR), pointer :: f_ptr(:)
    integer                      :: i

    ! Convert a C string (pointer to char) to a Fortran string

    l = INT(c_strlen(c_str))

    call C_F_POINTER(c_str, f_ptr, shape=[l])

    allocate(character(l) :: f_str)

    do i = 1, l
       f_str(i:i) = f_ptr(i)
    end do

    ! Finish

    return

  end function c_f_string

end module core_c
