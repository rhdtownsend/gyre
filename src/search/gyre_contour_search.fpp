! Module   : gyre_contour_search
! Purpose  : mode searching (complex, contour)
!
! Copyright 2013-2020 Rich Townsend & The GYRE Team
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

module gyre_contour_search

  ! Uses

  use core_kinds
  use core_memory
  use core_parallel
  use core_order

  use gyre_bvp
  use gyre_contour_map
  use gyre_discrim_func
  use gyre_ext
  use gyre_num_par
  use gyre_prox_search
  use gyre_state
  use gyre_status
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  integer, save :: j_m = 0

  ! Access specifiers

  private

  public :: contour_search

contains

  subroutine contour_search (bp, omega_re, omega_im, omega_min, omega_max, i, nm_p, process_mode)

    class(c_bvp_t), target, intent(inout) :: bp
    real(WP), intent(in)                  :: omega_re(:)
    real(WP), intent(in)                  :: omega_im(:)
    real(WP), intent(in)                  :: omega_min
    real(WP), intent(in)                  :: omega_max
    integer, intent(in)                   :: i
    type(num_par_t), intent(in)           :: nm_p
    interface
       subroutine process_mode (md, n_iter, chi)
         use core_kinds
         use gyre_ext
         use gyre_mode
         type(mode_t), intent(in)  :: md
         integer, intent(in)       :: n_iter
         type(r_ext_t), intent(in) :: chi
       end subroutine process_mode
    end interface

    type(c_state_t)          :: st
    type(c_discrim_func_t)   :: df
    complex(WP), allocatable :: omega_in_a(:)
    complex(WP), allocatable :: omega_in_b(:)
    integer, allocatable     :: j_in(:)

    ! Set up the discriminant function

     st = c_state_t()
     df = c_discrim_func_t(bp, st, omega_min, omega_max)

    ! Find contour intersections

    call find_isects_(df, omega_re, omega_im, i, nm_p, omega_in_a, omega_in_b, j_in)

    ! Search for modes

    call prox_search(bp, omega_in_a, omega_in_b, j_in, omega_min, omega_max, nm_p, process_mode)

    ! Finish

    return

  end subroutine contour_search

  !****

  subroutine find_isects_ (df, omega_re, omega_im, i, nm_p, omega_in_a, omega_in_b, j_in)

    type(c_discrim_func_t), intent(inout) :: df
    real(WP), intent(in)                  :: omega_re(:)
    real(WP), intent(in)                  :: omega_im(:)
    integer, intent(in)                   :: i
    type(num_par_t), intent(in)           :: nm_p
    complex(WP), allocatable, intent(out) :: omega_in_a(:)
    complex(WP), allocatable, intent(out) :: omega_in_b(:)
    integer, allocatable, intent(out)     :: j_in(:)

    integer                           :: n_omega_re
    integer                           :: n_omega_im
    integer, allocatable              :: k_part(:)
    integer                           :: n_percent
    integer                           :: c_beg
    integer                           :: c_end
    integer                           :: c_rate
    integer                           :: k
    integer                           :: j_ri(2)
    complex(WP)                       :: omega
    type(c_ext_t)                     :: discrim
    integer                           :: i_percent
    integer                           :: status
    complex(WP), allocatable          :: discrim_map_f(:,:)
    integer, allocatable              :: discrim_map_e(:,:)
    $if ($MPI)
    integer                           :: p
    $endif
    type(c_ext_t), allocatable        :: discrim_map(:,:)
    type(contour_map_t)               :: cm
    integer                           :: n_in
    integer                           :: d_in
    
    ! Calculate the discriminant on the omega grid

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Evaluating discriminant'
100    format(A)
    endif

    n_omega_re = SIZE(omega_re)
    n_omega_im = SIZE(omega_im)

    allocate(discrim_map_f(n_omega_re,n_omega_im))
    allocate(discrim_map_e(n_omega_re,n_omega_im))

    allocate(k_part(MPI_SIZE+1))

    call partition_tasks(n_omega_re*n_omega_im, 1, k_part)

    n_percent = 0

    call SYSTEM_CLOCK(c_beg, c_rate)

    discrim_loop: do k = k_part(MPI_RANK+1), k_part(MPI_RANK+2)-1

       j_ri = index_nd(k, [n_omega_re,n_omega_im])

       omega = CMPLX(omega_re(j_ri(1)), omega_im(j_ri(2)), KIND=WP)

       call df%eval(c_ext_t(omega), discrim, status)
       $ASSERT(status == STATUS_OK,Invalid status)

       discrim_map_f(j_ri(1),j_ri(2)) = FRACTION(discrim)
       discrim_map_e(j_ri(1),j_ri(2)) = EXPONENT(discrim)

       if (check_log_level('DEBUG')) then
          i_percent = FLOOR(100._WP*REAL(k-k_part(MPI_RANK+1))/REAL(k_part(MPI_RANK+2)-k_part(MPI_RANK+1)-1))
          if (i_percent > n_percent) then
             write(OUTPUT_UNIT, 110) 'Percent complete: ', i_percent
110          format(A,1X,I0)
             n_percent = i_percent
          end if
       endif

    end do discrim_loop

    $if ($MPI)

    do p = 1,MPI_SIZE
       call bcast_seq(discrim_map_f, k_part(p), k_part(p+1)-1, p-1)
       call bcast_seq(discrim_map_e, k_part(p), k_part(p+1)-1, p-1)
    end do

    $endif

    discrim_map = scale(c_ext_t(discrim_map_f), discrim_map_e)

    call SYSTEM_CLOCK(c_end)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'Time elapsed :', REAL(c_end-c_beg, WP)/c_rate, 's'
120    format(2X,A,F10.3,1X,A/)
    endif

    ! Create the contour map

    d_in = 128
    n_in = 0

    allocate(omega_in_a(d_in))
    allocate(omega_in_b(d_in))
    allocate(j_in(d_in))

    cm = contour_map_t(omega_re, omega_im, discrim_map, process_isect_)

    omega_in_a = omega_in_a(:n_in)
    omega_in_b = omega_in_b(:n_in)
    j_in = j_in(:n_in)

    ! Dump the map

    $if ($HDF5)
    call cm%dump(nm_p, i)
    $endif

    ! Finish

    return

  contains

    subroutine process_isect_ (omega_a_re, omega_b_re, omega_a_im, omega_b_im)

      complex(WP), intent(in) :: omega_a_re
      complex(WP), intent(in) :: omega_b_re
      complex(WP), intent(in) :: omega_a_im
      complex(WP), intent(in) :: omega_b_im

      ! Process the intersection

      n_in = n_in + 1

      if (n_in > d_in) then
         d_in = 2*d_in
         call reallocate(omega_in_a, [d_in])
         call reallocate(omega_in_b, [d_in])
         call reallocate(j_in, [d_in])
      endif

       j_m = j_m + 1

      omega_in_a(n_in) = omega_init_(omega_a_re, omega_b_re, 'im')
      omega_in_b(n_in) = omega_init_(omega_a_im, omega_b_im, 're')

      j_in(n_in) = j_m
    
      ! Finish

      return

    end subroutine process_isect_

    !****

    function omega_init_ (omega_a, omega_b, part) result (omega_init)

      complex(WP), intent(in)  :: omega_a
      complex(WP), intent(in)  :: omega_b
      character(*), intent(in) :: part
      complex(WP)              :: omega_init

      type(c_ext_t) :: discrim_a
      type(c_ext_t) :: discrim_b
      integer       :: status
      real(WP)      :: w

      call df%eval(c_ext_t(omega_a), discrim_a, status)
      $ASSERT(status == STATUS_OK,Invalid status)

      call df%eval(c_ext_t(omega_b), discrim_b, status)
      $ASSERT(status == STATUS_OK,Invalid status)
      
      ! Look for the point on the segment where the real/imaginary
      ! part of the discriminant changes sign

      select case (part)
      case ('re')
         w = real(-real_part(discrim_a)/(real_part(discrim_b) - real_part(discrim_a)))
      case ('im')
         w = real(-imag_part(discrim_a)/(imag_part(discrim_b) - imag_part(discrim_a)))
      case default
         $ABORT(Invalid part)
      end select

      omega_init = (1._WP-w)*omega_a + w*omega_b

      ! Finish

      return

    end function omega_init_

  end subroutine find_isects_

end module gyre_contour_search
