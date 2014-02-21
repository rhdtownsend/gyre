! Module   : gyre_famdl_file
! Purpose  : read formatted AMDL files
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

module gyre_famdl_file

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_coeffs
  use gyre_coeffs_evol
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_famdl_file

  ! Procedures

contains

  subroutine read_famdl_file (file, deriv_type, data_format, cf, x)

    character(LEN=*), intent(in)                 :: file
    character(LEN=*), intent(in)                 :: deriv_type
    character(LEN=*), intent(in)                 :: data_format
    type(coeffs_evol_t), intent(out)             :: cf
    real(WP), allocatable, intent(out), optional :: x(:)

    character(LEN=:), allocatable :: data_format_
    integer                       :: unit
    integer                       :: nmod
    integer                       :: n
    integer                       :: ivar
    real(WP), allocatable         :: glob(:)
    real(WP), allocatable         :: var(:,:)
    real(WP)                      :: M_star
    real(WP)                      :: R_star
    real(WP)                      :: L_star
    real(WP), allocatable         :: x_(:)
    real(WP), allocatable         :: c_1(:)
    real(WP), allocatable         :: V_g(:)
    real(WP), allocatable         :: Gamma_1(:) 
    real(WP), allocatable         :: As(:) 
    real(WP), allocatable         :: U(:)
    logical                       :: add_center

    if(data_format /= '') then
       data_format_ = data_format
    else
       data_format_ = '(1P4E20.13)'
    endif

    ! Read the model from the formatted AMDL-format file

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from FAMDL file', TRIM(file)
100    format(A,1X,A)
    endif

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header

    read(unit, 110) nmod, n, ivar
110 format(3I10)

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'Initial points :', n
120    format(2X,A,1X,I0)
    endif

    ! Read the data

    allocate(glob(8))
    allocate(var(ivar+1,n))

    read(unit, data_format_) glob, var

    close(unit)

    M_star = glob(1)
    R_star = glob(2)
    L_star = 0._WP

    x_ = var(1,:)

    c_1 = 1._WP/var(2,:)
    V_g = var(3,:)
    Gamma_1 = var(4,:)
    As = var(5,:)
    U = var(6,:)

    add_center = x_(1) /= 0._WP

    if(add_center .AND. check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Adding central point'
    endif

    ! Initialize the coeffs

    cf = coeffs_evol_t(M_star, R_star, L_star, x_, V_g*Gamma_1, As, U, c_1, Gamma_1, &
                       deriv_type, add_center)

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

  end subroutine read_famdl_file

end module gyre_famdl_file
