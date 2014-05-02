! Module   : gyre_losc_file
! Purpose  : read LOSC files
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

module gyre_losc_file

  ! Uses

  use core_kinds
  use core_order

  use gyre_constants
  use gyre_model
  use gyre_evol_model
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_losc_model

  ! Procedures

contains

  subroutine read_losc_model (file, deriv_type, regularize, ml, x)

    character(*), intent(in)                     :: file
    character(*), intent(in)                     :: deriv_type
    logical, intent(in)                          :: regularize
    type(evol_model_t), intent(out)              :: ml
    real(WP), allocatable, optional, intent(out) :: x(:)

    integer                   :: unit
    character(256)            :: line
    integer                   :: n
    real(WP)                  :: glob(3)
    real(WP), allocatable     :: var(:,:)
    integer                   :: i
    integer                   :: k
    integer, allocatable      :: ind(:)
    real(WP)                  :: M_star
    real(WP)                  :: R_star
    real(WP)                  :: L_star
    real(WP), allocatable     :: x_(:)
    real(WP), allocatable     :: c_1(:)
    real(WP), allocatable     :: V(:)
    real(WP), allocatable     :: Gamma_1(:) 
    real(WP), allocatable     :: As(:) 
    real(WP), allocatable     :: U(:)
    logical                   :: add_center

    ! Read data from the LOSC-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from LOSC file', TRIM(file)
100    format(A,1X,A)
    endif
          
    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header

    header_loop : do

       read(unit, *) line

       if (line == '%%beginoscdata') exit header_loop

    end do header_loop

    read(unit, *) glob

    read(unit, *)
    read(unit, *)
    read(unit, *) n

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Initial points :', n
110    format(3X,A,1X,I0)
    endif

    ! Read the data

    allocate(var(6,n))

    read_loop : do i = 1,n
       read(unit, *) k, var(:,i)
    end do read_loop

    close(unit)

    ind = unique_indices(var(1,:))

    if (SIZE(ind) < n) then

       if(check_log_level('WARN')) then
          write(OUTPUT_UNIT, 120) 'WARNING: Duplicate r-point(s) found, using innermost value(s)'
120       format('!!',1X,A)
       endif

       n = SIZE(var, 2)

    endif
       
    var = var(:,ind)

    M_star = glob(2)
    R_star = glob(1)
    L_star = 0._WP

    G_gravity = glob(3)

    x_ = var(1,:)/R_star

    c_1 = M_star/(R_star**3*var(2,:))
    V = var(4,:)*G_GRAVITY*var(2,:)*var(1,:)**2/var(3,:)
    U = 4._WP*PI*var(4,:)/var(2,:)
    Gamma_1 = var(5,:)
    As = -var(6,:)*var(1,:)**2

    add_center = x_(1) /= 0._WP

    if (check_log_level('INFO')) then
       if (add_center) write(OUTPUT_UNIT, 130) 'Adding central point'
130    format(3X,A)
    endif

    if (regularize) write(OUTPUT_UNIT, 130) 'Warning: Cannot regularize this model'

    ! Initialize the model

    ml = evol_model_t(M_star, R_star, L_star, x_, V, As, U, c_1, Gamma_1, &
                      deriv_type, add_center=add_center)

    ! Set up the grid

    if(PRESENT(x)) then
       if(add_center) then
          x = [0._WP,x_]
       else
          x = x_
       endif
    endif

    ! Finish

    return

  end subroutine read_losc_model

end module gyre_losc_file
