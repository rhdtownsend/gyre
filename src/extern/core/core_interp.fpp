! Module   : core_interp
! Purpose  : interpolant abstract type

$include 'core.inc'

module core_interp

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: interp_t
     integer :: n
   contains
     private
     procedure(f_n_), deferred     :: x_n_
     generic, public               :: x => x_n_
     procedure(f_1_), deferred     :: y_1_
     procedure(f_v_), deferred     :: y_v_
     procedure(f_n_), deferred     :: y_n_
     generic, public               :: y => y_1_, y_v_, y_n_
     procedure(f_1_), deferred     :: dy_1_
     procedure(f_v_), deferred     :: dy_v_
     procedure(f_n_), deferred     :: dy_n_
     generic, public               :: dy => dy_1_, dy_v_, dy_n_
     procedure(f_1_1_), deferred   :: iy_1_1_
     procedure(f_1_v_), deferred   :: iy_1_v_
     procedure(f_v_1_), deferred   :: iy_v_1_
     procedure(f_v_v_), deferred   :: iy_v_v_
     generic, public               :: iy => iy_1_1_, iy_1_v_, iy_v_1_, iy_v_v_
     procedure(attribs_), deferred :: attribs_
     generic, public               :: attribs => attribs_
  end type interp_t

  ! Interfaces

  abstract interface

     function f_1_ (this, x) result (f)
       use core_kinds
       import interp_t
       class(interp_t), intent(in) :: this
       real(WP), intent(in)        :: x
       real(WP)                    :: f
     end function f_1_

     function f_v_ (this, x) result (f)
       use core_kinds
       import interp_t
       class(interp_t), intent(in) :: this
       real(WP), intent(in)        :: x(:)
       real(WP)                    :: f(SIZE(x))
     end function f_v_

     function f_n_ (this) result (f)
       use core_kinds
       import interp_t
       class(interp_t), intent(in) :: this
       real(WP)                    :: f(this%n)
     end function f_n_

     function f_1_1_ (this, x_a, x_b) result (f)
       use core_kinds
       import interp_t
       class(interp_t), intent(in) :: this
       real(WP), intent(in)        :: x_a
       real(WP), intent(in)        :: x_b
       real(WP)                    :: f
     end function f_1_1_

     function f_1_v_ (this, x_a, x_b) result (f)
       use core_kinds
       import interp_t
       class(interp_t), intent(in) :: this
       real(WP), intent(in)        :: x_a
       real(WP), intent(in)        :: x_b(:)
       real(WP)                    :: f(SIZE(x_b))
     end function f_1_v_

     function f_v_1_ (this, x_a, x_b) result (f)
       use core_kinds
       import interp_t
       class(interp_t), intent(in) :: this
       real(WP), intent(in)        :: x_a(:)
       real(WP), intent(in)        :: x_b
       real(WP)                    :: f(SIZE(x_a))
     end function f_v_1_

     function f_v_v_ (this, x_a, x_b) result (f)
       use core_kinds
       import interp_t
       class(interp_t), intent(in) :: this
       real(WP), intent(in)        :: x_a(:)
       real(WP), intent(in)        :: x_b(:)
       real(WP)                    :: f(SIZE(x_a))
     end function f_v_v_

     subroutine attribs_ (this, x_min, x_max, n)
       use core_kinds
       import interp_t
       class(interp_t), intent(in)     :: this
       real(WP), optional, intent(out) :: x_min
       real(WP), optional, intent(out) :: x_max
       integer, optional, intent(out)  :: n
     end subroutine attribs_

  end interface

  ! Access specifiers

  private

  public :: interp_t

end module core_interp
