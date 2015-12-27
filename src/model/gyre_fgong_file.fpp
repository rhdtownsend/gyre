! Module   : gyre_fgong_file
! Purpose  : read FGONG files
!
! Copyright 2015 Rich Townsend
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
  use gyre_model_par
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

    type(model_par_t), intent(in)   :: ml_p
    type(evol_model_t), intent(out) :: ml

    integer                   :: unit
    integer                   :: n
    integer                   :: iconst
    integer                   :: ivar
    integer                   :: ivers
    character(:), allocatable :: data_format
    real(WP), allocatable     :: glob(:)
    real(WP), allocatable     :: var(:,:)
    integer                   :: i
    real(WP)                  :: M_star
    real(WP)                  :: R_star
    real(WP)                  :: L_star
    real(WP), allocatable     :: r(:)
    real(WP), allocatable     :: m(:)
    real(WP), allocatable     :: p(:)
    real(WP), allocatable     :: rho(:) 
    real(WP), allocatable     :: T(:) 
    real(WP), allocatable     :: N2(:)
    real(WP), allocatable     :: Gamma_1(:)
    real(WP), allocatable     :: nabla_ad(:)
    real(WP), allocatable     :: delta(:)
    real(WP), allocatable     :: x(:)
    real(WP), allocatable     :: V_2(:)
    real(WP), allocatable     :: As(:)
    real(WP), allocatable     :: U(:)
    real(WP), allocatable     :: c_1(:)
    real(WP), allocatable     :: Omega_rot(:)
    logical                   :: has_center

    ! Open the FGONG-format file

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from FGONG file', TRIM(ml_p%file)
100    format(A,1X,A)
    endif

    open(NEWUNIT=unit, FILE=ml_p%file, STATUS='OLD')

    ! Read the header

    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)

    read(unit, *) n, iconst, ivar, ivers

     if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Initial points :', n
       write(OUTPUT_UNIT, 110) 'File version   :', ivers
110    format(3X,A,1X,I0)
    endif

    ! Read the data

    if (ml_p%data_format /= '') then
       data_format = ml_p%data_format
    else
       if (ivers < 1000) then
          data_format = '(1P5E16.9)'
       else
          $ABORT(Cannot handle ivers > 1000)
       endif
    endif

    allocate(glob(iconst))
    allocate(var(ivar,n))

    read(unit, data_format) glob

    read_loop : do i = 1, n
       read(unit, data_format) var(:,i)
    end do read_loop

    close(unit)

    ! Extract structure data

    M_star = glob(1)
    R_star = glob(2)
    L_star = glob(3)

    r = var(1,:)/R_star

    m = EXP(var(2,:))*M_star
    T = var(3,:)
    p = var(4,:)
    rho = var(5,:)

    Gamma_1 = var(10,:)
    nabla_ad = var(11,:)
    delta = var(12,:)

    allocate(N2(n))

    where (r/R_star >= EPSILON(0._WP))
       N2 = G_GRAVITY*m*var(15,:)/r**3
    elsewhere
       N2 = 0._WP
    endwhere

    if (r(1)/R_star < EPSILON(0._WP)) r(1) = 0._WP
    if (m(1)/M_star < EPSILON(0._WP)) m(1) = 0._WP

    if (m(1) == 0._WP .AND. r(1) /= 0._WP) then
       r(1) = 0._WP
       write(OUTPUT_UNIT, 130) 'Forcing central r == 0'
130    format(3X,A)
    elseif(r(1) == 0._WP .AND. m(1) /= 0._WP) then
       m(1) = 0._WP
       write(OUTPUT_UNIT, 130) 'Forcing central m == 0'
    endif

    ! Calculate dimensionless structure data

    allocate(x(n))

    allocate(V_2(n))
    allocate(As(n))
    allocate(U(n))
    allocate(c_1(n))

    allocate(Omega_rot(n))

    x = r/R_star

    where (x /= 0._WP)
       V_2 = G_GRAVITY*m*rho/(p*r*x**2)
       As = r**3*N2/(G_GRAVITY*m)
       U = 4._WP*PI*rho*r**3/m
       c_1 = (r/R_star)**3/(m/M_star)
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*p(1))
       As = 0._WP
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
    end where

    if (ml_p%uniform_rot) then
       Omega_rot = ml_p%Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    else
       Omega_rot = 0._WP
    endif

    ! Initialize the model

    ml = evol_model_t(x, M_star, R_star, L_star, ml_p)

    call ml%set_V_2(V_2)
    call ml%set_As(As)
    call ml%set_U(U)
    call ml%set_c_1(c_1)

    call ml%set_Gamma_1(Gamma_1)
    call ml%set_delta(delta)
    call ml%set_nabla_ad(nabla_ad)

    call ml%set_Omega_rot(Omega_rot)

    ! Finish

    return

  end subroutine read_fgong_model

end module gyre_fgong_file
