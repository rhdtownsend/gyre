! Module   : gyre_ad_oscpar
! Purpose  : adiabatic oscillation parameters

$include 'core.inc'

module gyre_ad_oscpar

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ad_oscpar_t
     private
     real(WP), public              :: lambda_0
     integer, public               :: l
!    character(LEN=:), allocatable :: outer_bound_type
     character(LEN=256), public    :: outer_bound_type
   contains
     private
     procedure, public :: init
  end type ad_oscpar_t

 ! Access specifiers

  private

  public :: ad_oscpar_t

  ! Procedures

contains

  subroutine init (this, l, outer_bound_type)

    class(ad_oscpar_t), intent(out) :: this
    integer, intent(in)             :: l
    character(LEN=*), intent(in)    :: outer_bound_type

    ! Initialize the ad_oscpar

    if(l == 0) then
       this%lambda_0 = 0._WP
    else
       this%lambda_0 = l - 2._WP
    endif

    this%l = l

    this%outer_bound_type = outer_bound_type

    ! Finish

    return

  end subroutine init

end module gyre_ad_oscpar
