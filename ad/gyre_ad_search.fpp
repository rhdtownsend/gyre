! Module   : gyre_ad_search
! Purpose  : adiabatic mode searching
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

module gyre_ad_search

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_numpar
  use gyre_ad_bvp
  use gyre_ad_discfunc
  use gyre_mode
  use gyre_frontend
  use gyre_ext_arith

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: ad_scan_search

  ! Procedures

contains

  subroutine ad_scan_search (bp, omega, md)

    type(ad_bvp_t), target, intent(inout)  :: bp
    real(WP), intent(in)                   :: omega(:)
    type(mode_t), allocatable, intent(out) :: md(:)

    integer             :: n_omega
    integer             :: i_part(MPI_SIZE+1)
    integer             :: c_beg
    integer             :: c_end
    integer             :: c_rate
    integer             :: i
    type(ext_real_t)    :: discrim(SIZE(omega))
    $if($MPI)
    integer             :: recvcounts(MPI_SIZE)
    integer             :: displs(MPI_SIZE)
    $endif
    integer             :: n_brack
    integer             :: i_brack(SIZE(omega)-1)
    type(ad_discfunc_t) :: df
    integer             :: n_iter
    complex(WP)         :: omega_root
    $if($MPI)
    integer             :: p
    $endif

    ! Calculate the discriminant on the omega abscissa

    call write_header('Adiabatic Discriminant Scan', '=')

    n_omega = SIZE(omega)

    call partition_tasks(n_omega, 1, i_part)

    call SYSTEM_CLOCK(c_beg, c_rate)

    $if($MPI)
    call barrier()
    $endif

    discrim_loop : do i = i_part(MPI_RANK+1),i_part(MPI_RANK+2)-1

       discrim(i) = ext_real(bp%discrim(CMPLX(omega(i), KIND=WP)))

       write(OUTPUT_UNIT, 100) 'Eval:', omega(i), discrim(i)%f, discrim(i)%e
100    format(A,2X,E23.16,2X,F19.16,2X,I7)

    end do discrim_loop

    $if($MPI)
    call barrier()
    $endif

    call SYSTEM_CLOCK(c_end)
    if(MPI_RANK == 0) then
       write(OUTPUT_UNIT, 110) 'Completed ad scan; time elapsed:', REAL(c_end-c_beg, WP)/c_rate, 's'
110    format(/A,1X,F10.3,1X,A)
    endif

    ! Gather data

    $if($MPI)

    recvcounts = i_part(2:)-i_part(:MPI_SIZE)
    displs = i_part(:MPI_SIZE)-1

    call allgatherv(discrim, recvcounts, displs)

    $endif

    ! Scan for root brackets

    n_brack = 0

    bracket_loop : do i = 1, n_omega-1

       if(discrim(i)*discrim(i+1) <= ext_real(0._WP)) then

          ! Store the bracket location

          n_brack = n_brack + 1
          i_brack(n_brack) = i

       end if

    end do bracket_loop

    ! Process each bracket to find modes

    call write_header('Adiabatic Mode Finding', '=')

    call df%init(bp)

    call partition_tasks(n_brack, 1, i_part)

    allocate(md(n_brack))

    call SYSTEM_CLOCK(c_beg, c_rate)

    $if($MPI)
    call barrier()
    $endif

    root_loop : do i = i_part(MPI_RANK+1), i_part(MPI_RANK+2)-1

       ! Set the discriminant normalization, based on the mid-bracket
       ! frequency

       call bp%set_norm(CMPLX(0.5_WP*(omega(i_brack(i)) + omega(i_brack(i)+1)), KIND=WP))

       ! Find the root

       n_iter = bp%np%n_iter_max

       omega_root = df%root(omega(i_brack(i)), omega(i_brack(i)+1), 0._WP, n_iter=n_iter)
       $ASSERT(n_iter <= bp%np%n_iter_max,Too many iterations)

       ! Set up the mode

       call md(i)%init(bp, omega_root)

       ! Report

       write(OUTPUT_UNIT, '(A,3(2X,I6),3(2X,E23.16),2X,I4)') 'Mode :', md(i)%n_p-md(i)%n_g, md(i)%n_p, md(i)%n_g, &
            md(i)%omega, ABS(md(i)%discrim), n_iter

    end do root_loop

    $if($MPI)
    call barrier()
    $endif

    call SYSTEM_CLOCK(c_end)
    if(MPI_RANK == 0) then
       write(OUTPUT_UNIT, 110) 'Completed ad find; time elapsed:', REAL(c_end-c_beg, WP)/c_rate, 's'
    endif

    ! Broadcast data

    $if($MPI)

    do p = 0, MPI_SIZE-1
       do i = i_part(p+1), i_part(p+2)-1
          call md(i)%bcast(p)
       end do
    enddo

    $endif

    ! Finish

    return

  end subroutine ad_scan_search

end module gyre_ad_search
