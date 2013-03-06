! Module   : gyre_jacobian
! Purpose  : Jacobian evaluation (interface)

$include 'core.inc'

module gyre_jacobian

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: jacobian_t
     private
     integer, public :: n_e
   contains
     private
     procedure(eval_i), deferred, public :: eval
     procedure(eval_i), deferred, public :: eval_logx
  end type jacobian_t

  ! Interfaces

  abstract interface

     subroutine eval_i (this, omega, x, A)
       use core_kinds
       import jacobian_t
       class(jacobian_t), intent(in) :: this
       complex(WP), intent(in)       :: omega
       real(WP), intent(in)          :: x
       complex(WP), intent(out)      :: A(:,:)
     end subroutine eval_i

  end interface

  ! Access specifiers

  private

  public :: jacobian_t

end module gyre_jacobian
