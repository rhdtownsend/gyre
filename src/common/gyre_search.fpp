! Module   : gyre_search
! Purpose  : mode searching
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

module gyre_search

  ! Uses

  use core_kinds
  use core_constants
  use core_order
  use core_parallel

  use gyre_bvp
  use gyre_coeffs
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_scanpar
  use gyre_mode
  use gyre_grid
  use gyre_ext_arith
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: build_scan
  public :: scan_search
  public :: prox_search

contains

  subroutine build_scan (sp, cf, op, gp, x_in, omega)

    type(scanpar_t), intent(in)        :: sp(:)
    class(coeffs_t), intent(in)        :: cf
    type(oscpar_t), intent(in)         :: op
    type(gridpar_t), intent(in)        :: gp(:)
    real(WP), allocatable, intent(in)  :: x_in(:)
    real(WP), allocatable, intent(out) :: omega(:)

    real(WP) :: x_i
    real(WP) :: x_o
    integer  :: i
    real(WP) :: omega_min
    real(WP) :: omega_max
    integer  :: j

    ! Determine the grid range

    call grid_range(gp, cf, op, x_in, x_i, x_o)

    ! Loop through scanpars

    allocate(omega(0))

    sp_loop : do i = 1,SIZE(sp)

       ! Set up the frequency grid

       omega_min = sp(i)%freq_min/freq_scale(cf, op, x_o, sp(i)%freq_units)
       omega_max = sp(i)%freq_max/freq_scale(cf, op, x_o, sp(i)%freq_units)
       
       select case(sp(i)%grid_type)
       case('LINEAR')
          omega = [omega,(((sp(i)%n_freq-j)*omega_min + (j-1)*omega_max)/(sp(i)%n_freq-1), j=1,sp(i)%n_freq)]
       case('INVERSE')
          omega = [omega,((sp(i)%n_freq-1)/((sp(i)%n_freq-j)/omega_min + (j-1)/omega_max), j=1,sp(i)%n_freq)]
       case default
          $ABORT(Invalid grid_type)
       end select

    end do sp_loop

    ! Sort the frequencies

    omega = omega(sort_indices(omega))

    ! Finish

    return

  end subroutine build_scan

!****

  subroutine scan_search (bp, omega, md)

    class(bvp_t), target, intent(inout)    :: bp
    real(WP), intent(in)                   :: omega(:)
    type(mode_t), allocatable, intent(out) :: md(:)

    real(WP), allocatable         :: omega_a(:)
    real(WP), allocatable         :: omega_b(:)
    type(ext_real_t), allocatable :: discrim_a(:)
    type(ext_real_t), allocatable :: discrim_b(:)
    integer                       :: n_brack
    integer                       :: i_part(MPI_SIZE+1)
    integer                       :: c_beg
    integer                       :: c_end
    integer                       :: c_rate
    integer                       :: i
    integer                       :: n_p
    integer                       :: n_g
    integer                       :: n_pg
    $if($MPI)
    integer                       :: p
    $endif

    ! Find discriminant root brackets

    call find_brackets(bp, omega, omega_a, omega_b, discrim_a, discrim_b)

    ! Process each bracket to find modes

    if(check_log_level('INFO')) then

       write(OUTPUT_UNIT, 100) form_header('Scan Mode Search', '=')
100    format(A)

       write(OUTPUT_UNIT, 110) 'l', 'n_pg', 'n_p', 'n_g', 'Re(omega)', 'Im(omega)', 'chi', 'n_iter'
110    format(4(2X,A6),3(2X,A24),2X,A6)
       
    endif

    n_brack = SIZE(omega_a)

    call partition_tasks(n_brack, 1, i_part)

    allocate(md(n_brack))

    call SYSTEM_CLOCK(c_beg, c_rate)

    $if($MPI)
    call barrier()
    $endif

    mode_loop : do i = i_part(MPI_RANK+1), i_part(MPI_RANK+2)-1

       ! Find the mode

       md(i) = bp%mode(CMPLX([omega_a(i),omega_b(i)], KIND=WP), &
                       ext_complex([discrim_a(i),discrim_b(i)]), use_real=.TRUE.)

       ! Report

       call md(i)%classify(n_p, n_g, n_pg)

       if(check_log_level('INFO', MPI_RANK)) then
          write(OUTPUT_UNIT, 120) md(i)%op%l, n_pg, n_p, n_g, md(i)%omega, real(md(i)%chi), md(i)%n_iter
120       format(4(2X,I6),3(2X,E24.16),2X,I6)
       endif

    end do mode_loop

    $if($MPI)
    call barrier()
    $endif

    call SYSTEM_CLOCK(c_end)

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Completed mode search; time elapsed:', REAL(c_end-c_beg, WP)/c_rate, 's'
130    format(/A,1X,F10.3,1X,A)
    endif

    ! Broadcast data

    $if($MPI)

    do p = 0, MPI_SIZE-1
       do i = i_part(p+1), i_part(p+2)-1
          call bcast(md(i), p)
       end do
    enddo

    $endif

    ! Finish

    return

  end subroutine scan_search

!****

  subroutine prox_search (bp, md)

    class(bvp_t), target, intent(inout) :: bp
    type(mode_t), intent(inout)         :: md(:)

    integer     :: n_md
    integer     :: i_part(MPI_SIZE+1)
    integer     :: c_beg
    integer     :: c_end
    integer     :: c_rate
    integer     :: i
    complex(WP) :: omega_a
    complex(WP) :: omega_b
    integer     :: n_p
    integer     :: n_g
    integer     :: n_pg
    $if($MPI)
    integer     :: p
    $endif

    ! Process each mode to find a proximate mode

    if(check_log_level('INFO')) then

       write(OUTPUT_UNIT, 100) form_header('Proximity Mode Search', '=')
100    format(A)

       write(OUTPUT_UNIT, 110) 'l', 'n_pg', 'n_p', 'n_g', 'Re(omega)', 'Im(omega)', 'chi', 'n_iter'
110    format(4(2X,A6),3(2X,A23),2X,A4)
       
    endif

    n_md = SIZE(md)

    call partition_tasks(n_md, 1, i_part)

    call SYSTEM_CLOCK(c_beg, c_rate)

    $if($MPI)
    call barrier()
    $endif

    mode_loop : do i = i_part(MPI_RANK+1), i_part(MPI_RANK+2)-1

       ! Find the mode

       omega_a = md(i)%omega*CMPLX(1._WP,  SQRT(EPSILON(0._WP)), WP)
       omega_b = md(i)%omega*CMPLX(1._WP, -SQRT(EPSILON(0._WP)), WP)

       md(i) = bp%mode([omega_a,omega_b])

       ! Report

       call md(i)%classify(n_p, n_g, n_pg)

       if(check_log_level('INFO', MPI_RANK)) then
          write(OUTPUT_UNIT, 120) md(i)%op%l, n_pg, n_p, n_g, md(i)%omega, real(md(i)%chi), md(i)%n_iter
120       format(4(2X,I6),3(2X,E24.16),2X,I4)
       endif

    end do mode_loop

    $if($MPI)
    call barrier()
    $endif

    call SYSTEM_CLOCK(c_end)

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Completed mode search; time elapsed:', REAL(c_end-c_beg, WP)/c_rate, 's'
130    format(/A,1X,F10.3,1X,A)
    endif

    ! Broadcast data

    $if($MPI)

    do p = 0, MPI_SIZE-1
       do i = i_part(p+1), i_part(p+2)-1
          call bcast(md(i), p)
       end do
    enddo

    $endif

    ! Finish

    return

  end subroutine prox_search

!****

  subroutine find_brackets (bp, omega, omega_a, omega_b, discrim_a, discrim_b)

    class(bvp_t), target, intent(inout)        :: bp
    real(WP), intent(in)                       :: omega(:)
    real(WP), allocatable, intent(out)         :: omega_a(:)
    real(WP), allocatable, intent(out)         :: omega_b(:)
    type(ext_real_t), allocatable, intent(out) :: discrim_a(:)
    type(ext_real_t), allocatable, intent(out) :: discrim_b(:)

    integer          :: n_omega
    integer          :: i_part(MPI_SIZE+1)
    integer          :: c_beg
    integer          :: c_end
    integer          :: c_rate
    integer          :: i
    type(ext_real_t) :: discrim(SIZE(omega))
    $if($MPI)
    integer          :: recvcounts(MPI_SIZE)
    integer          :: displs(MPI_SIZE)
    $endif
    integer          :: n_brack
    integer          :: i_brack(SIZE(omega))

    ! Calculate the discriminant on the omega abscissa

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) form_header('Discriminant Scan', '=')
100    format(A)
    endif

    n_omega = SIZE(omega)

    call partition_tasks(n_omega, 1, i_part)

    call SYSTEM_CLOCK(c_beg, c_rate)

    $if($MPI)
    call barrier()
    $endif

    discrim_loop : do i = i_part(MPI_RANK+1),i_part(MPI_RANK+2)-1

       discrim(i) = ext_real(bp%discrim(CMPLX(omega(i), KIND=WP)))

       if(check_log_level('DEBUG', MPI_RANK)) then
          write(OUTPUT_UNIT, 110) 'Eval:', omega(i), fraction(discrim(i)), exponent(discrim(i))
110       format(A,2X,E24.16,2X,F19.16,2X,I7)
       endif

    end do discrim_loop

    $if($MPI)

    recvcounts = i_part(2:)-i_part(:MPI_SIZE)
    displs = i_part(:MPI_SIZE)-1

    call allgatherv(discrim, recvcounts, displs)

    $endif

    $if($MPI)
    call barrier()
    $endif

    call SYSTEM_CLOCK(c_end)

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'Completed scan; time elapsed:', REAL(c_end-c_beg, WP)/c_rate, 's'
120    format(/A,1X,F10.3,1X,A)
    endif

    ! Find root brackets

    n_brack = 0

    bracket_loop : do i = 1, n_omega-1

       if(discrim(i)*discrim(i+1) <= ext_real(0._WP)) then
          n_brack = n_brack + 1
          i_brack(n_brack) = i
       end if

    end do bracket_loop

    ! Set up the bracket frequencies

    omega_a = omega(i_brack(:n_brack))
    omega_b = omega(i_brack(:n_brack)+1)

    discrim_a = discrim(i_brack(:n_brack))
    discrim_b = discrim(i_brack(:n_brack)+1)

    ! Finish

    return

  end subroutine find_brackets

end module gyre_search
