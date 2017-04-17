! Module   : gyre_losc_file
! Purpose  : read LOSC files
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

module gyre_losc_file

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

  public :: read_losc_model

  ! Procedures

contains

  subroutine read_losc_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    integer                     :: unit
    character(256)              :: line
    integer                     :: n
    real(WP)                    :: glob(3)
    real(WP), allocatable       :: var(:,:)
    integer                     :: i
    integer                     :: k
    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP), allocatable       :: x(:)
    real(WP), allocatable       :: V_2(:)
    real(WP), allocatable       :: As(:) 
    real(WP), allocatable       :: U(:)
    real(WP), allocatable       :: c_1(:)
    real(WP), allocatable       :: Gamma_1(:) 
    real(WP), allocatable       :: Omega_rot(:)
    type(evol_model_t), pointer :: em

    ! Open the LOSC-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from LOSC file'
100    format(A)
       write(OUTPUT_UNIT, 110) 'File name', TRIM(ml_p%file)
110    format(3X,A,1X,A)
    endif
          
    open(NEWUNIT=unit, FILE=ml_p%file, STATUS='OLD')

    ! Read the header

    header_loop : do

       read(unit, *) line

       if (line == '%%beginoscdata') exit header_loop

    end do header_loop

    read(unit, *) glob

    read(unit, *)
    read(unit, *)
    read(unit, *) n

    ! Read the data

    allocate(var(6,n))

    read_loop : do i = 1,n
       read(unit, *) k, var(:,i)
    end do read_loop

    close(unit)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'Read', n, 'points'
120    format(3X,A,1X,I0,1X,A)
    endif
    
    ! Extract structure data

    M_star = glob(2)
    R_star = glob(1)

    if (glob(3) /= G_GRAVITY) then
       $WARN(Gravitational constant in LOSC file does not match value specified in &constants)
    endif

    x = var(1,:)/R_star

    allocate(V_2(n))

    where (x /= 0._WP)
       V_2 = G_GRAVITY*var(4,:)*var(2,:)*R_star**2/var(3,:)
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*var(4,:)**2*R_star**2/(3._WP*var(3,:))
    endwhere
       
    As = -var(6,:)*var(1,:)**2
    U = 4._WP*PI*var(4,:)/var(2,:)
    c_1 = M_star/(R_star**3*var(2,:))

    Gamma_1 = var(5,:)

    ! Snap grid points

    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), x)

    ! Calculate dimensionless structure data

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

  end subroutine read_losc_model

end module gyre_losc_file
