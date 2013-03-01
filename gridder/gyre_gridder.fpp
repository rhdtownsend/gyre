! Module   : gyre_gridder
! Purpose  : grid construction (interface)

$include 'core.inc'

module gyre_gridder

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: gridder_t
   contains
     private
     procedure(build_i), deferred, public :: build
  end type gridder_t

  ! Interfaces

  abstract interface

     subroutine build_i (this, omega, x)
       use core_kinds
       import gridder_t
       class(gridder_t), intent(in)       :: this
       complex(WP), intent(in)            :: omega
       real(WP), intent(out), allocatable :: x(:)
     end subroutine build_i

  end interface

  ! Access specifiers

  private

  public :: gridder_t

end module gyre_gridder
