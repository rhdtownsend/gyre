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

  use gyre_numpar
  use gyre_nad_bvp
  use gyre_nad_discfunc
  use gyre_eigfunc
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

  subroutine nad_prox_search (bp, ad_ef, nad_ef)

    type(nad_bvp_t), target, intent(inout)    :: bp
    type(eigfunc_t), intent(in)               :: ad_ef(:)
    type(eigfunc_t), allocatable, intent(out) :: nad_ef(:)

    integer                 :: n_ef
    integer                 :: i_part(MPI_SIZE+1)
    integer                 :: c_beg
    integer                 :: c_end
    integer                 :: c_rate
    integer                 :: i
    type(nad_discfunc_t)    :: df
    type(numpar_t), pointer :: np
    integer                 :: n_iter
    complex(WP)             :: omega_root
    complex(WP)             :: discrim_root
    integer                 :: n_p
    integer                 :: n_g
    $if($MPI)
    integer                 :: p
    $endif

    ! Process each adiabatic root to find non-adiabatic root

    call write_header('Non-Adiabatic Mode Finding', '=')

    call df%init(bp)

    n_ef = SIZE(ad_ef)

    call partition_tasks(n_ef, 1, i_part)

    allocate(nad_ef(n_ef))

    np => bp%get_np()

    call SYSTEM_CLOCK(c_beg, c_rate)

    $if($MPI)
    call barrier()
    $endif

    root_loop : do i = i_part(MPI_RANK+1), i_part(MPI_RANK+2)-1

       ! Set the discriminant normalization, based on the adiabatic
       ! frequency

       call bp%set_x_ad(ad_ef(i)%omega)
       call bp%set_norm(ad_ef(i)%omega)

       ! Find the root

       n_iter = np%n_iter_max

       omega_root = df%root(ad_ef(i)%omega*CMPLX(1._WP, SQRT(EPSILON(0._WP)), WP), &
                            ad_ef(i)%omega*CMPLX(1._WP, -SQRT(EPSILON(0._WP)), WP), &
                            0._WP, n_iter=n_iter)
       $ASSERT(n_iter <= np%n_iter_max,Too many iterations)

       discrim_root = df%eval(omega_root)

       ! Set up the eigfunction

       nad_ef(i) = bp%eigfunc(omega_root)

       ! Report

       call nad_ef(i)%classify(n_p, n_g)

       write(OUTPUT_UNIT, 100) 'Mode :', n_p-n_g, n_p, n_g, omega_root, ABS(discrim_root), n_iter
100    format(A,3(2X,I6),3(2X,E23.16),2X,I4)

    end do root_loop

    $if($MPI)
    call barrier()
    $endif

    call SYSTEM_CLOCK(c_end)
    if(MPI_RANK == 0) then
       write(OUTPUT_UNIT, 110) 'Completed nad find; time elapsed:', REAL(c_end-c_beg, WP)/c_rate, 's'
110    format(/A,1X,F10.3,1X,A)
    endif

    ! Broadcast data

    $if($MPI)

    do p = 0, MPI_SIZE-1
       do i = i_part(p+1), i_part(p+2)-1
          call bcast(nad_ef(i), p)
       end do
    enddo

    $endif

    ! Finish

    return

  end subroutine nad_prox_search

end module gyre_nad_search
