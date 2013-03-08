! Module   : gyre_nad_oscpar
! Purpose  : nonadiabatic oscillation parameters

$include 'core.inc'

module gyre_nad_oscpar

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: nad_oscpar_t
     private
     real(WP), public              :: lambda_0
     real(WP), public              :: ad_thresh
     integer, public               :: l
!    character(LEN=:), allocatable :: outer_bound_type
     character(LEN=256), public    :: outer_bound_type
     logical, public               :: force_ad
   contains
     private
     procedure, public :: init
  end type nad_oscpar_t

 ! Access specifiers

  private

  public :: nad_oscpar_t

  ! Procedures

contains

  subroutine init (this, l, outer_bound_type, force_ad, ad_thresh)

    class(nad_oscpar_t), intent(out) :: this
    integer, intent(in)              :: l
    character(LEN=*), intent(in)     :: outer_bound_type
    logical, intent(in)              :: force_ad
    real(WP), intent(in)             :: ad_thresh

    ! Initialize the nad_oscpar

    if(l == 0) then
       this%lambda_0 = 0._WP
    else
       this%lambda_0 = l - 2._WP
    endif

    this%l = l

    this%outer_bound_type = outer_bound_type

    this%force_ad = force_ad
    this%ad_thresh = ad_thresh

    ! Finish

    return

  end subroutine init

end module gyre_nad_oscpar
