! Program  : gyre_writer
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
     private
     procedure(init_i), deferred, public  :: init
     procedure(final_i), deferred, public :: final
     procedure(write_i_0_i), deferred     :: write_i_0
     procedure(write_i_1_i), deferred     :: write_i_1
     procedure(write_r_0_i), deferred     :: write_r_0
     procedure(write_r_1_i), deferred     :: write_r_1
     procedure(write_c_0_i), deferred     :: write_c_0
     procedure(write_c_1_i), deferred     :: write_c_1
     procedure(write_a_0_i), deferred     :: write_a_0
     procedure(write_a_1_i), deferred     :: write_a_1
     generic, public                      :: write => write_i_0, write_i_1, write_r_0, write_r_1, &
                                                      write_c_0, write_c_1, write_a_0, write_a_1
  end type writer_t

  ! Interfaces

  abstract interface

     subroutine init_i (this, file_name)
       import writer_t
       class(writer_t), intent(out) :: this
       character(LEN=*), intent(in) :: file_name
     end subroutine init_i

     subroutine final_i (this)
       import writer_t
       class(writer_t), intent(inout) :: this
     end subroutine final_i

     $define $WRITE_I $sub

     $local $INFIX $1
     $local $DATA_TYPE $2
     $local $DATA_RANK $3

     subroutine write_${INFIX}_${DATA_RANK}_i (this, name, data)
       use core_kinds
       import writer_t
       class(writer_t), intent(inout) :: this
       character(LEN=*), intent(in)   :: name
       $DATA_TYPE, intent(in)         :: data$ARRAY_SPEC($DATA_RANK)
     end subroutine write_${INFIX}_${DATA_RANK}_i

     $endsub

     $WRITE_I(i,integer,0)
     $WRITE_I(i,integer,1)

     $WRITE_I(r,real(WP),0)
     $WRITE_I(r,real(WP),1)

     $WRITE_I(c,complex(WP),0)
     $WRITE_I(c,complex(WP),1)

     $WRITE_I(a,character(LEN=*),0)
     $WRITE_I(a,character(LEN=*),1)

  end interface

  ! Access specifiers

  private

  public :: writer_t

end module gyre_writer
