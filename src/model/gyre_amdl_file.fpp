! Module   : gyre_amdl_file
! Purpose  : read unformatted AMDL files
!
! Copyright 2013-2014 Rich Townsend
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

module gyre_amdl_file

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

  public :: read_amdl_model

  ! Procedures

contains

  subroutine read_amdl_model (file, deriv_type, add_center, ml, x)

    character(*), intent(in)                     :: file
    character(*), intent(in)                     :: deriv_type
    logical, intent(in)                          :: add_center
    type(evol_model_t), intent(out)              :: ml
    real(WP), allocatable, intent(out), optional :: x(:)

    integer                   :: unit
    integer                   :: nmod
    integer                   :: n
    integer                   :: ivar
    real(WP)                  :: glob(8)
    integer                   :: idata8
    real(WP), allocatable     :: var(:,:)
    integer, allocatable      :: ind(:)
    real(WP)                  :: M_star
    real(WP)                  :: R_star
    real(WP)                  :: L_star
    real(WP), allocatable     :: x_(:)
    real(WP), allocatable     :: c_1(:)
    real(WP), allocatable     :: V_g(:)
    real(WP), allocatable     :: Gamma_1(:) 
    real(WP), allocatable     :: As(:) 
    real(WP), allocatable     :: U(:)
    logical                   :: has_center

    ! Read the model from the unformatted AMDL-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from AMDL file', TRIM(file)
100    format(A,1X,A)
    endif

    open(NEWUNIT=unit, FILE=file, STATUS='OLD', FORM='UNFORMATTED')

    ! Read the header

    read(unit) nmod, n, glob

    idata8 = int(glob(8)+0.1_WP)

    if (idata8 >= 100) then
       ivar = 8
    elseif (idata8 >= 10) then
       ivar = 6
    else
       ivar = 5
    endif
       
    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'Initial points :', n
120    format(3X,A,1X,I0)
    endif

    ! Read the data

    rewind(unit)

    allocate(var(ivar+1,n))

    read(unit) nmod, n, glob, var

    close(unit)

    ind = unique_indices(var(1,:))

    if (SIZE(ind) < n) then

       if(check_log_level('WARN')) then
          write(OUTPUT_UNIT, 130) 'WARNING: Duplicate x-point(s) found, using innermost value(s)'
130       format('!!',1X,A)
       endif

       n = SIZE(var, 2)

    endif
       
    var = var(:,ind)

    M_star = glob(1)
    R_star = glob(2)
    L_star = 0._WP

    x_ = var(1,:)

    c_1 = 1._WP/var(2,:)
    V_g = var(3,:)
    Gamma_1 = var(4,:)
    As = var(5,:)
    U = var(6,:)

    has_center = x_(1) == 0._WP

    if (check_log_level('INFO')) then
       if (add_center) then
          if (has_center) then
             write(OUTPUT_UNIT, 140) 'No need to add central point'
140          format(3X,A)
          else
             write(OUTPUT_UNIT, 140) 'Adding central point'
          endif
       endif
    endif

    ! Initialize the model

    ml = evol_model_t(M_star, R_star, L_star, x_, V_g*Gamma_1, As, U, c_1, Gamma_1, &
                      deriv_type, add_center=add_center .AND. .NOT. has_center)

    ! Set up the grid

    if (PRESENT(x)) then
       if (add_center .AND. .NOT. has_center) then
          x = [0._WP,x_]
       else
          x = x_
       endif
    endif

    ! Finish

    return

  end subroutine read_amdl_model

end module gyre_amdl_file
