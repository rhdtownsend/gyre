! Module   : gyre_mode
! Purpose  : mode classification

$include 'core.inc'

module gyre_mode

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_bvp
  use gyre_ext_arith

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: mode_t
!     private
     real(WP), allocatable    :: x(:)
     complex(WP), allocatable :: y(:,:)
     complex(WP), public      :: omega
     complex(WP), public      :: discrim
     integer                  :: n
     integer, public          :: n_p
     integer, public          :: n_g
   contains
     private
     procedure, public :: init
     $if($MPI)
     procedure, public :: bcast => bcast_md
     $endif
  end type mode_t

  ! Access specifiers

  private

  public :: mode_t

  ! Procedures

contains

  subroutine init (this, bp, omega)

    class(mode_t), intent(out)  :: this
    class(bvp_t), intent(inout) :: bp
    complex(WP), intent(in)     :: omega

    ! Initialize the mode

    this%omega = omega
    this%discrim = cmplx(bp%discrim(omega))

    call bp%recon(omega, this%x, this%y)

    call classify(this%x, REAL(this%y), this%n_p, this%n_g)

    this%n = SIZE(this%x)

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_md (this, root_rank)

    class(mode_t), intent(inout) :: this
    integer, intent(in)          :: root_rank

    ! Broadcast the mode

    call bcast(this%x, root_rank, alloc=.TRUE.)
    call bcast(this%y, root_rank, alloc=.TRUE.)

    call bcast(this%omega, root_rank)

    call bcast(this%n, root_rank)

    call bcast(this%n_p, root_rank)
    call bcast(this%n_g, root_rank)

  end subroutine bcast_md

  $endif

!****

  subroutine classify (x, y, n_p, n_g)

    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: y(:,:)
    integer, intent(out) :: n_p
    integer, intent(out) :: n_g

    logical  :: inner_ext
    integer  :: j
    real(WP) :: y_2_cross

    $ASSERT(SIZE(y, 1) >= 2,Dimension mismatch)
    $ASSERT(SIZE(y, 2) == SIZE(x),Dimension mismatch)

    ! Classify the mode using the Cowling-Scuflaire scheme

    n_p = 0
    n_g = 0

    inner_ext = ABS(y(1,1)) > ABS(y(1,2))

    x_loop : do j = 2,SIZE(x)-1

       ! If the innermost extremum in y_1 hasn't yet been reached,
       ! skip

       if(.NOT. inner_ext) then
          inner_ext = ABS(y(1,j)) > ABS(y(1,j-1)) .AND. ABS(y(1,j)) > ABS(y(1,j+1))
          cycle x_loop
       endif

       ! Look for a node in y_1

       if(y(1,j) >= 0._WP .AND. y(1,j+1) < 0._WP) then

          y_2_cross = y(2,j) - y(1,j)*(y(2,j+1) - y(2,j))/(y(1,j+1) - y(1,j))

          if(y_2_cross >= 0._WP) then
             n_p = n_p + 1
          else
             n_g = n_g + 1
          endif

       elseif(y(1,j) <= 0._WP .AND. y(1,j+1) > 0._WP) then

         y_2_cross = y(2,j) - y(1,j)*(y(2,j+1) - y(2,j))/(y(1,j+1) - y(1,j))

          if(y_2_cross <= 0._WP) then
             n_p = n_p + 1
          else
             n_g = n_g + 1
          endif

       endif

    end do x_loop

    ! Finish

    return

  end subroutine classify

end module gyre_mode
