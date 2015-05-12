! Module   : gyre_fgong_file
! Purpose  : read FGONG files
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

module gyre_fgong_file

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

  public :: read_fgong_model

  ! Procedures

contains

  subroutine read_fgong_model (file, deriv_type, data_format, add_center, ml, x)

    character(*), intent(in)                     :: file
    character(*), intent(in)                     :: deriv_type
    character(*), intent(in)                     :: data_format
    logical, intent(in)                          :: add_center
    type(evol_model_t), intent(out)              :: ml
    real(WP), allocatable, intent(out), optional :: x(:)

    character(:), allocatable :: data_format_
    integer                   :: unit
    integer                   :: n
    integer                   :: iconst
    integer                   :: ivar
    integer                   :: ivers
    real(WP), allocatable     :: glob(:)
    real(WP), allocatable     :: var(:,:)
    integer                   :: i
    integer, allocatable      :: ind(:)
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
    logical                   :: has_center

    if(data_format /= '') then
       data_format_ = data_format
    else
       data_format_ = '(1P5E16.9)'
    endif

    ! Read data from the FGONG-format file

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from FGONG file', TRIM(file)
100    format(A,1X,A)
    endif

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

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

    allocate(glob(iconst))
    allocate(var(ivar,n))

    read(unit, data_format_) glob

    read_loop : do i = 1,n
       read(unit, data_format_) var(:,i)
    end do read_loop

    close(unit)

    ind = unique_indices(var(1,:))

    if (SIZE(ind) < n) then

       if(check_log_level('WARN')) then
          write(OUTPUT_UNIT, 120) 'WARNING: Duplicate x-point(s) found, using innermost value(s)'
120       format('!!',1X,A)
       endif

       n = SIZE(ind)

    endif

    var = var(:,ind)

    M_star = glob(1)
    R_star = glob(2)
    L_star = glob(3)

    r = var(1,:)
    m = EXP(var(2,:))*M_star
    T = var(3,:)
    p = var(4,:)
    rho = var(5,:)
    Gamma_1 = var(10,:)
    nabla_ad = var(11,:)
    delta = var(12,:)

    allocate(N2(n))

    where(r/R_star >= EPSILON(0._WP))
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
    elseif (r(1) == 0._WP .AND. m(1) /= 0._WP) then
       m(1) = 0._WP
       write(OUTPUT_UNIT, 130) 'Forcing central m == 0'
    endif

    has_center = r(1) == 0._WP .AND. m(1) == 0._WP

    if (check_log_level('INFO')) then
       if (add_center) then
          if (has_center) then
             write(OUTPUT_UNIT, 130) 'No need to add central point'
          else
             write(OUTPUT_UNIT, 130) 'Adding central point'
          endif
       endif
    endif

    ! Initialize the model

    ml = evol_model_t(M_star, R_star, L_star, r, m, p, rho, T, &
                      N2, Gamma_1, nabla_ad, delta, SPREAD(0._WP, DIM=1, NCOPIES=n), &
                      deriv_type, add_center=add_center .AND. .NOT. has_center)

    ! Set up the grid

    if (PRESENT(x)) then
       if (add_center .AND. .NOT. has_center) then
          x = [0._WP,r/R_star]
       else
          x = r/R_star
       endif
    endif

    ! Finish

    return

  end subroutine read_fgong_model

end module gyre_fgong_file
