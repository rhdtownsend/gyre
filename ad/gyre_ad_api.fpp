! Module   : gyre_ad_api
! Purpose  : Adiabatic API

$include 'core.inc'

module gyre_ad_api

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_ad_bvp
  use gyre_mode
  use gyre_util
  use gyre_ext_arith

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: find_ad_modes

  ! Procedures

contains

  subroutine find_ad_modes (bp, omega, n_iter_max, md)

    type(ad_bvp_t), target, intent(inout)  :: bp
    real(WP), intent(in)                   :: omega(:)
    integer, intent(in)                    :: n_iter_max
    type(mode_t), allocatable, intent(out) :: md(:)

    integer                  :: n_omega
    integer                  :: i_part(MPI_SIZE+1)
    integer                  :: i
    type(ext_real_t)         :: discrim(SIZE(omega))
    $if($MPI)
    integer                  :: recvcounts(MPI_SIZE)
    integer                  :: displs(MPI_SIZE)
    $endif
    integer                  :: n_brack
    integer                  :: i_brack(SIZE(omega)-1)
    integer                  :: n_iter
    complex(WP)              :: omega_root
    real(WP), allocatable    :: x(:)
    complex(WP), allocatable :: y(:,:)
    $if($MPI)
    integer            :: p
    $endif

    ! Calculate the discriminant on the omega abscissa

    call write_header('Adiabatic Discriminant Scan', '=')

    n_omega = SIZE(omega)

    call partition_tasks(n_omega, 1, i_part)

    discrim_loop : do i = i_part(MPI_RANK+1),i_part(MPI_RANK+2)-1

       discrim(i) = bp%discrim(omega(i))

       write(OUTPUT_UNIT, *) 'Eval:',omega(i),discrim(i)%f,discrim(i)%e

    end do discrim_loop

    ! Gather data

    $if($MPI)

    recvcounts = i_part(2:)-i_part(:MPI_SIZE)
    displs = i_part(:MPI_SIZE)-1

    call allgatherv(discrim, recvcounts, displs)

    $endif

    ! Search for root brackets

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

    call partition_tasks(n_brack, 1, i_part)

    allocate(md(n_brack))

    root_loop : do i = i_part(MPI_RANK+1), i_part(MPI_RANK+2)-1

       ! Set up the normalizing exponent

       call bp%set_norm(CMPLX(0.5_WP*(omega(i_brack(i)) + omega(i_brack(i)+1)), KIND=WP))

       ! Find the root

       n_iter = n_iter_max

       omega_root = bp%root(omega(i_brack(i)), omega(i_brack(i)+1), 0._WP, n_iter=n_iter)
       $ASSERT(n_iter <= n_iter_max,Too many iterations)

       ! Set up the mode

       call bp%recon(omega_root, x, y)
       call md(i)%init(omega_root, bp%eval(omega_root), x, y)

       ! Report

       write(OUTPUT_UNIT, '(A,3(2X,I6),3(2X,E23.16),2X,I4)') 'Mode :', md(i)%n_p-md(i)%n_g, md(i)%n_p, md(i)%n_g, &
            md(i)%omega, ABS(md(i)%discrim), n_iter

    end do root_loop

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

  end subroutine find_ad_modes

end module gyre_ad_api
