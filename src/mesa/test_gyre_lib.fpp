! Program  : test_gyre_lib
! Purpose  : test gyre_lib for memory leaks and other issues
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

program test_gyre_lib

  ! Uses

  use core_kinds
  use gyre_constants

  use gyre_lib

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  real(WP) :: rpar(1)
  integer  :: ipar(1)
  integer  :: i
  
  ! Initialize

  call gyre_init('gyre.in')

  ! Read the model

  call gyre_read_model('bcep.mesa')

  ! Repeatedy find modes

  do i = 1,100
     print *,'Iteration:',i
     call gyre_get_modes(0, user_sub, ipar, rpar)
     call gyre_get_modes(1, user_sub, ipar, rpar)
     call gyre_get_modes(2, user_sub, ipar, rpar)
  end do

  ! Finish

contains

  subroutine user_sub (md, ipar, rpar, retcode)

    use core_kinds
    use gyre_mode
    type(mode_t), intent(in) :: md
    integer, intent(inout)   :: ipar(:)
    real(WP), intent(inout)  :: rpar(:)
    integer, intent(out)     :: retcode
    
    ! Print out mode data

    print *, md%n_pg, md%omega

    ! Finish

    retcode = 0

    return

  end subroutine user_sub

end program test_gyre_lib


  

  
