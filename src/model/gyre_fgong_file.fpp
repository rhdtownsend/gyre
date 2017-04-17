! Module   : gyre_fgong_file
! Purpose  : read FGONG files
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

module gyre_fgong_file

  ! Uses

  use core_kinds

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

  public :: read_fgong_model

  ! Procedures

contains

  subroutine read_fgong_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    integer                     :: unit
    integer                     :: n
    integer                     :: iconst
    integer                     :: ivar
    integer                     :: ivers
    character(:), allocatable   :: data_format
    real(WP), allocatable       :: glob(:)
    real(WP), allocatable       :: var(:,:)
    integer                     :: i
    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: L_star
    real(WP), allocatable       :: x(:)
    real(WP), allocatable       :: m(:)
    real(WP), allocatable       :: P(:)
    real(WP), allocatable       :: rho(:) 
    real(WP), allocatable       :: T(:) 
    real(WP), allocatable       :: Gamma_1(:)
    real(WP), allocatable       :: nabla_ad(:)
    real(WP), allocatable       :: delta(:)
    real(WP), allocatable       :: V_2(:)
    real(WP), allocatable       :: As(:)
    real(WP), allocatable       :: U(:)
    real(WP), allocatable       :: c_1(:)
    real(WP), allocatable       :: Omega_rot(:)
    type(evol_model_t), pointer :: em

    ! Open the FGONG-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from FGONG file'
100    format(A)
       write(OUTPUT_UNIT, 110) 'File name', TRIM(ml_p%file)
110    format(3X,A,1X,A)
    endif

    open(NEWUNIT=unit, FILE=ml_p%file, STATUS='OLD')

    ! Read the header

    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)

    read(unit, *) n, iconst, ivar, ivers

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'File version', ivers
120    format(3X,A,1X,I0)
    endif

    ! Read the data

    if (ml_p%data_format /= '') then
       data_format = ml_p%data_format
    else
       if (ivers < 1000) then
          data_format = '(1P5E16.9)'
       else
          data_format = '(1P,5(X,E26.18E3))'
       endif
    endif

    allocate(glob(iconst))
    allocate(var(ivar,n))

    read(unit, data_format) glob

    read_loop : do i = 1, n
       read(unit, data_format) var(:,i)
    end do read_loop

    close(unit)

    var = var(:,n:1:-1)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Read', n, 'points'
130    format(3X,A,1X,I0,1X,A)
    endif
    
    ! Extract structure data

    M_star = glob(1)
    R_star = glob(2)
    L_star = glob(3)

    x = var(1,:)/R_star

    $if ($GFORTRAN_PR_69185)
    allocate(m(n))
    $endif

    m = EXP(var(2,:))
    T = var(3,:)
    P = var(4,:)
    rho = var(5,:)

    Gamma_1 = var(10,:)
    nabla_ad = var(11,:)
    delta = var(12,:)

    As = var(15,:)

    ! Snap grid points

    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), x, m)

    ! Calculate dimensionless structure data

    allocate(V_2(n))
    allocate(U(n))
    allocate(c_1(n))

    where (x /= 0._WP)
       V_2 = G_GRAVITY*(m*M_star)*rho/(P*x**3*R_star)
       U = 4._WP*PI*rho*(x*R_star)**3/(m*M_star)
       c_1 = x**3/m
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*P(1))
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
    end where

    allocate(Omega_rot(n))

    if (ml_p%uniform_rot) then
       Omega_rot = uniform_Omega_rot(ml_p, M_star, R_star)
    else
       Omega_rot = 0._WP
    endif

    ! Initialize the evol_model_t

    allocate(em, SOURCE=evol_model_t(x, M_star, R_star, L_star, ml_p))

    call em%define(I_V_2, V_2)
    call em%define(I_AS, As)
    call em%define(I_U, U)
    call em%define(I_C_1, c_1)

    call em%define(I_GAMMA_1, Gamma_1)
    call em%define(I_DELTA, delta)
    call em%define(I_NABLA_AD, nabla_ad)

    call em%define(I_OMEGA_ROT, Omega_rot)

    ! Return a pointer

    ml => em

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine read_fgong_model

end module gyre_fgong_file
