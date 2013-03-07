! Module   : gyre_ad_params
! Purpose  : adiabatic oscillation parameters

$include 'core.inc'

module gyre_ad_params

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ad_params_t
     private
     real(WP), public :: lambda_0
     integer, public  :: l
   contains
     private
     procedure, public :: init
     $if($MPI)
     procedure, public :: bcast => bcast_ap
     $endif
  end type ad_params_t

 ! Access specifiers

  private

  public :: ad_params_t

  ! Procedures

contains

  subroutine init (this, l)

    class(ad_params_t), intent(out) :: this
    integer, intent(in)             :: l

    ! Initialize the ad_params

    if(l == 0) then
       this%lambda_0 = 0._WP
    else
       this%lambda_0 = l - 2._WP
    endif

    this%l = l

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_ap (this, root_rank)

    class(ad_params_t), intent(inout) :: this
    integer, intent(in)               :: root_rank

    ! Broadcast the ad_params

    call bcast(this%lambda_0, root_rank)
    call bcast(this%l, root_rank)

    ! Finish

    return

  end subroutine bcast_ap

  $endif

end module gyre_ad_params
