! Module   : gyre_nad_search
! Purpose  : nonadiabatic mode search
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

module gyre_nad_search

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_nad_bvp
  use gyre_nad_discfunc
  use gyre_mode
  use gyre_frontend
  use gyre_ext_arith

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: nad_prox_search

  ! Procedures

contains

  subroutine nad_prox_search (bp, ad_md, n_iter_max, nad_md)

    type(nad_bvp_t), target, intent(inout)  :: bp
    type(mode_t), intent(in)                :: ad_md(:)
    integer, intent(in)                     :: n_iter_max
    type(mode_t), allocatable, intent(out)  :: nad_md(:)

    integer              :: n_md
    integer              :: i_part(MPI_SIZE+1)
    integer              :: c_beg
    integer              :: c_end
    integer              :: c_rate
    integer              :: i
    type(nad_discfunc_t) :: df
    integer              :: n_iter
    complex(WP)          :: omega_root
    $if($MPI)
    integer              :: p
    $endif

    ! Process each adiabatic mode to find non-adiabatic modes

    call write_header('Non-Adiabatic Mode Finding', '=')

    call df%init(bp)

    n_md = SIZE(ad_md)

    call partition_tasks(n_md, 1, i_part)

    allocate(nad_md(n_md))

    call SYSTEM_CLOCK(c_beg, c_rate)

    root_loop : do i = i_part(MPI_RANK+1), i_part(MPI_RANK+2)-1

       ! Set the discriminant normalization, based on the adiabatic
       ! frequency

       call bp%set_norm(ad_md(i)%omega)

       ! Find the root

       n_iter = n_iter_max

       omega_root = df%root(ad_md(i)%omega*CMPLX(1._WP, SQRT(EPSILON(0._WP)), WP), &
                            ad_md(i)%omega*CMPLX(1._WP, -SQRT(EPSILON(0._WP)), WP), &
                            0._WP, n_iter=n_iter)
       $ASSERT(n_iter <= n_iter_max,Too many iterations)

       ! Set up the mode

       call nad_md(i)%init(bp, omega_root)

       ! Report

       write(OUTPUT_UNIT, '(A,3(2X,I6),3(2X,E23.16),2X,I4)') 'Mode :', nad_md(i)%n_p-nad_md(i)%n_g, nad_md(i)%n_p, nad_md(i)%n_g, &
            nad_md(i)%omega, ABS(nad_md(i)%discrim), n_iter

    end do root_loop

    call SYSTEM_CLOCK(c_end)
    if(MPI_RANK == 0) then
       write(OUTPUT_UNIT, 100) 'Completed nad find; time elapsed:', REAL(c_end-c_beg, WP)/c_rate, 's'
100    format(/A,1X,F10.3,1X,A)
    endif

    ! Broadcast data

    $if($MPI)

    do p = 0, MPI_SIZE-1
       do i = i_part(p+1), i_part(p+2)-1
          call nad_md(i)%bcast(p)
       end do
    enddo

    $endif

    ! Finish

    return

  end subroutine nad_prox_search

end module gyre_nad_search
