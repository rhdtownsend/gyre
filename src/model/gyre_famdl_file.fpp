! Module   : gyre_famdl_file
! Purpose  : read FAMDL files
!
! Copyright 2013-2017 Rich Townsend
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
  use core_order

  use gyre_constants
  use gyre_evol_model
  use gyre_model
  use gyre_model_par
  use gyre_model_util
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_famdl_model

  ! Procedures

contains

  subroutine read_famdl_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    character(:), allocatable   :: data_format
    integer                     :: unit
    integer                     :: nmod
    integer                     :: n
    integer                     :: ivar
    real(WP), allocatable       :: glob(:)
    real(WP), allocatable       :: var(:,:)
    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: P_c
    real(WP)                    :: rho_c
    real(WP), allocatable       :: x(:)
    real(WP), allocatable       :: c_1(:)
    real(WP), allocatable       :: V_g(:)
    real(WP), allocatable       :: Gamma_1(:) 
    real(WP), allocatable       :: As(:) 
    real(WP), allocatable       :: U(:)
    real(WP), allocatable       :: V_2(:)
    real(WP), allocatable       :: Omega_rot(:)
    type(evol_model_t), pointer :: em

    ! Open the AMDL-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from FAMDL file'
100    format(A)
       write(OUTPUT_UNIT, 110) 'File name', TRIM(ml_p%file)
110    format(3X,A,1X,A)
    endif

    open(NEWUNIT=unit, FILE=ml_p%file, STATUS='OLD')

    ! Read the header

    read(unit, 120) nmod, n, ivar
120 format(3I10)

    ! Read the data

    if (ml_p%data_format /= '') then
       data_format = ml_p%data_format
    else
       data_format = '(1P4E20.13)'
    endif

    allocate(glob(8))
    allocate(var(ivar+1,n))

    read(unit, data_format) glob, var

    close(unit)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Read', n, 'points'
130    format(3X,A,1X,I0,1X,A)
    endif

    ! Extract structure data

    M_star = glob(1)
    R_star = glob(2)
    P_c = glob(3)
    rho_c = glob(4)

    x = var(1,:)

    c_1 = 1._WP/var(2,:)
    V_g = var(3,:)
    Gamma_1 = var(4,:)
    As = var(5,:)
    U = var(6,:)

    ! Snap grid points

    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), x)
  
    ! Calculate dimensionless structure data

    allocate(V_2(n))

    where (x /= 0._WP)
       V_2 = V_g*Gamma_1
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho_c**2*R_star**2/(3._WP*P_c)
    end where

    allocate(Omega_rot(n))

    if (ml_p%uniform_rot) then
       Omega_rot = uniform_Omega_rot(ml_p, M_star, R_star)
    else
       Omega_rot = 0._WP
    endif

    ! Initialize the evol_model_t

    allocate(em, SOURCE=evol_model_t(x, M_star, R_star, 0._WP, ml_p))

    call em%define(I_V_2, V_2)
    call em%define(I_AS, As)
    call em%define(I_U, U)
    call em%define(I_C_1, c_1)

    call em%define(I_GAMMA_1, Gamma_1)

    call em%define(I_OMEGA_ROT, Omega_rot)

    ! Return a pointer

    ml => em

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine read_famdl_model

end module gyre_famdl_file
