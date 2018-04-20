! Module   : gyre_writer
! Purpose  : write data (interface)
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

module gyre_writer

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: writer_t
   contains
     procedure(final), deferred     :: final
     procedure(write_i_0), deferred :: write_i_0
     procedure(write_i_1), deferred :: write_i_1
     procedure(write_r_0), deferred :: write_r_0
     procedure(write_r_1), deferred :: write_r_1
     procedure(write_c_0), deferred :: write_c_0
     procedure(write_c_1), deferred :: write_c_1
     procedure(write_a_0), deferred :: write_a_0
     procedure(write_a_1), deferred :: write_a_1
     generic, public                :: write => write_i_0, write_i_1, write_r_0, write_r_1, &
                                                write_c_0, write_c_1, write_a_0, write_a_1
  end type writer_t

  ! Interfaces

  abstract interface

     subroutine final (this)
       import writer_t
       class(writer_t), intent(inout) :: this
     end subroutine final

     $define $WRITE_ $sub

     $local $INFIX $1
     $local $DATA_TYPE $2
     $local $DATA_RANK $3

     subroutine write_${INFIX}_${DATA_RANK} (this, name, data)
       use core_kinds
       import writer_t
       class(writer_t), intent(inout) :: this
       character(*), intent(in)       :: name
       $DATA_TYPE, intent(in)         :: data$ARRAY_SPEC($DATA_RANK)
     end subroutine write_${INFIX}_${DATA_RANK}

     $endsub

     $WRITE_(i,integer,0)
     $WRITE_(i,integer,1)

     $WRITE_(r,real(WP),0)
     $WRITE_(r,real(WP),1)

     $WRITE_(c,complex(WP),0)
     $WRITE_(c,complex(WP),1)

     $WRITE_(a,character(*),0)
     $WRITE_(a,character(*),1)

  end interface

  ! Access specifiers

  private

  public :: writer_t

end module gyre_writer
