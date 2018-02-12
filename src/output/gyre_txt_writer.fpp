! Module   : gyre_txt_writer
! Purpose  : write txt (ASCII) data
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

module gyre_txt_writer

  ! Uses

  use core_kinds
  use core_memory

  use gyre_writer
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: FIELD_LEN = 25

  character(*), parameter :: I_FORMAT = '(I25)'
  character(*), parameter :: R_FORMAT = '(E25.16E3)'
  character(*), parameter :: A_FORMAT = '(A25)'

  ! Derived-type definitions

  type, extends (writer_t) :: txt_writer_t
     private
     character(FIELD_LEN), allocatable :: s_data(:)
     character(FIELD_LEN), allocatable :: v_data(:,:)
     character(FIELD_LEN), allocatable :: s_names(:)
     character(FIELD_LEN), allocatable :: v_names(:)
     integer                           :: unit
     integer                           :: n_s
     integer                           :: n_v
     integer                           :: n
   contains
     private
     procedure, public :: final
     procedure, public :: write_i_0
     procedure, public :: write_i_1
     procedure, public :: write_r_0
     procedure, public :: write_r_1
     procedure, public :: write_c_0
     procedure, public :: write_c_1
     procedure, public :: write_a_0
     procedure, public :: write_a_1
  end type txt_writer_t

  ! Interfaces

  interface txt_writer_t
     module procedure txt_writer_t_
  end interface txt_writer_t

  ! Access specifiers

  private

  public :: txt_writer_t

  ! Procedures

contains

  function txt_writer_t_ (file_name, label) result (wr)

    character(*), intent(in) :: file_name
    character(*), intent(in) :: label
    type(txt_writer_t)       :: wr

    ! Construct the txt_writer_t

    open(NEWUNIT=wr%unit, FILE=file_name, STATUS='REPLACE', FORM='FORMATTED', ACCESS='STREAM')
    
    write(wr%unit, 100) label
100 format(A)

    wr%n_s = 0
    wr%n_v = 0

    wr%n = 0

    ! Finish

    return

  end function txt_writer_t_

!****

  subroutine final (this)

    class(txt_writer_t), intent(inout) :: this

    character(:), allocatable :: fmt
    integer                   :: i
    integer                   :: k

    ! Finalize the txt_writer

    fmt = join_fmts([I_FORMAT],[this%n_s])

    write(this%unit, fmt) [(i,i=1,this%n_s)]

    fmt = join_fmts([A_FORMAT],[this%n_s])

    write(this%unit, fmt) [(this%s_names(i),i=1,this%n_s)]
    write(this%unit, fmt) [(this%s_data(i)//'|',i=1,this%n_s)]
                                                   
    fmt = join_fmts([I_FORMAT],[this%n_v])

    write(this%unit, fmt) [(i,i=1,this%n_v)]

    fmt = join_fmts([A_FORMAT],[this%n_v])

    write(this%unit, fmt) [(this%v_names(i),i=1,this%n_v)]

    do k = 1, this%n
       write(this%unit, fmt) [(this%v_data(k,i)//'|',i=1,this%n_v)]
    end do

    close(this%unit)

    ! Finish

    return

  end subroutine final

!****
  
  $define $WRITE $sub

  $local $INFIX $1
  $local $DATA_TYPE $2
  $local $FORMAT $3

  subroutine write_${INFIX}_0 (this, name, data)

    class(txt_writer_t), intent(inout) :: this
    character(*), intent(in)           :: name
    $DATA_TYPE, intent(in)             :: data

    integer :: d

    ! If necessary, allocate/reallocate arrays

    if(this%n_s == 0) then
       allocate(this%s_data(16))
       allocate(this%s_names(16))
    endif

    this%n_s = this%n_s + 1

    d = SIZE(this%s_names)

    if(this%n_s > d) then
       call reallocate(this%s_names, [2*d])
       call reallocate(this%s_data, [2*d])
    endif

    ! Store the data

    this%s_names(this%n_s) = rjust(name, FIELD_LEN)

    write(this%s_data(this%n_s), $FORMAT) data

    ! Finish

    return

  end subroutine write_${INFIX}_0

!****

  subroutine write_${INFIX}_1 (this, name, data)

    class(txt_writer_t), intent(inout) :: this
    character(*), intent(in)           :: name
    $DATA_TYPE, intent(in)             :: data(:)

    integer :: d
    integer :: k

    ! If necessary, allocate/reallocate arrays

    if(this%n_v == 0) then
       this%n = SIZE(data)
       allocate(this%v_data(this%n,16))
       allocate(this%v_names(16))
    else
       $CHECK_BOUNDS(SIZE(data),this%n)
    endif

    this%n_v = this%n_v + 1

    d = SIZE(this%v_names)

    if(this%n_v > d) then
       call reallocate(this%v_names, [2*d])
       call reallocate(this%v_data, [this%n,2*d])
    endif

    ! Store the data

    this%v_names(this%n_v) = rjust(name, FIELD_LEN)

    do k = 1, this%n
       write(this%v_data(k,this%n_v), $FORMAT) data(k)
    end do

    ! Finish

    return

  end subroutine write_${INFIX}_1

  $endsub

  $WRITE(i,integer,I_FORMAT)
  $WRITE(r,real(WP),R_FORMAT)
  $WRITE(a,character(*),A_FORMAT)

!****
  
  subroutine write_c_0 (this, name, data)

    class(txt_writer_t), intent(inout) :: this
    character(*), intent(in)           :: name
    complex(WP), intent(in)            :: data

    integer :: d

    ! If necessary, allocate/reallocate arrays

    if(this%n_s == 0) then
       allocate(this%s_data(16))
       allocate(this%s_names(16))
    endif

    this%n_s = this%n_s + 2

    d = SIZE(this%s_names)

    if(this%n_s > d) then
       call reallocate(this%s_names, [2*d])
       call reallocate(this%s_data, [2*d])
    endif

    ! Store the data

    this%s_names(this%n_s-1) = rjust('Re('//TRIM(name)//')', FIELD_LEN)
    this%s_names(this%n_s  ) = rjust('Im('//TRIM(name)//')', FIELD_LEN)

    write(this%s_data(this%n_s-1), R_FORMAT) REAL(data)
    write(this%s_data(this%n_s  ), R_FORMAT) AIMAG(data)

    ! Finish

    return

  end subroutine write_c_0

!****

  subroutine write_c_1 (this, name, data)

    class(txt_writer_t), intent(inout) :: this
    character(*), intent(in)           :: name
    complex(WP), intent(in)            :: data(:)

    integer :: d
    integer :: k

    ! If necessary, allocate/reallocate arrays

    if(this%n_v == 0) then
       this%n = SIZE(data)
       allocate(this%v_data(this%n,16))
       allocate(this%v_names(16))
    else
       $CHECK_BOUNDS(SIZE(data),this%n)
    endif

    this%n_v = this%n_v + 2

    d = SIZE(this%v_names)

    if(this%n_v > d) then
       call reallocate(this%v_names, [2*d])
       call reallocate(this%v_data, [this%n,2*d])
    endif

    ! Store the data

    this%v_names(this%n_v-1) = rjust('Re('//TRIM(name)//')', FIELD_LEN)
    this%v_names(this%n_v  ) = rjust('Im('//TRIM(name)//')', FIELD_LEN)

    do k = 1, this%n
       write(this%v_data(k,this%n_v-1), R_FORMAT) REAL(data(k))
       write(this%v_data(k,this%n_v  ), R_FORMAT) AIMAG(data(k))
    end do

    ! Finish

    return

  end subroutine write_c_1

end module gyre_txt_writer
