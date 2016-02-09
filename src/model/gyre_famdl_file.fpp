! Module   : gyre_famdl_file
! Purpose  : read FAMDL files
!
! Copyright 2013-2015 Rich Townsend
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
       write(OUTPUT_UNIT, 100) 'Reading from FAMDL file', TRIM(ml_p%file)
100    format(A,1X,A)
    endif

    open(NEWUNIT=unit, FILE=ml_p%file, STATUS='OLD')

    ! Read the header

    read(unit, 110) nmod, n, ivar
110 format(3I10)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'Grid points  :', n
120    format(3X,A,1X,I0)
    endif

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

    ! Extract structure data

    M_star = glob(1)
    R_star = glob(2)
    P_c = glob(3)
    rho_c = glob(4)

<<<<<<< local
    p_c = glob(3)
    rho_c = glob(4)

    x_ = var(1,:)
=======
    x = var(1,:)
>>>>>>> other

    c_1 = 1._WP/var(2,:)
    V_g = var(3,:)
    Gamma_1 = var(4,:)
    As = var(5,:)
    U = var(6,:)

    ! Snap grid points

<<<<<<< local
    allocate(V_2(n))

    where (x_ /= 0._WP)
       V_2 = V_g*Gamma_1/x_**2
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho_c**2*R_star**2/(3._WP*p_c)
    end where

    if (check_log_level('INFO')) then
       if (add_center) then
          if (has_center) then
             write(OUTPUT_UNIT, 140) 'No need to add central point'
140          format(3X,A)
          else
             write(OUTPUT_UNIT, 140) 'Adding central point'
          endif
       endif
=======
    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), x)
  
    ! Calculate dimensionless structure data

    allocate(V_2(n))

    allocate(Omega_rot(n))

    where (x /= 0._WP)
       V_2 = V_g*Gamma_1
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho_c**2*R_star**2/(3._WP*P_c)
    end where

    if (ml_p%uniform_rot) then
       Omega_rot = ml_p%Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    else
       Omega_rot = 0._WP
>>>>>>> other
    endif

    ! Initialize the evol_model_t

<<<<<<< local
    ml = evol_model_t(M_star, R_star, L_star, x_, V_2, As, U, c_1, Gamma_1, &
                      deriv_type, add_center=add_center .AND. .NOT. has_center)
=======
    allocate(em, SOURCE=evol_model_t(x, M_star, R_star, 0._WP, ml_p))
>>>>>>> other

    call em%set_V_2(V_2)
    call em%set_As(As)
    call em%set_U(U)
    call em%set_c_1(c_1)

    call em%set_Gamma_1(Gamma_1)

    call em%set_Omega_rot(Omega_rot)

    ! Return a pointer

    ml => em

    ! Finish

    return

  end subroutine read_famdl_model

end module gyre_famdl_file
