! Module   : gyre_therm_coeffs
! Purpose  : thermal structure coefficients (interface)

$include 'core.inc'

module gyre_therm_coeffs

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure(get_1_i), deferred :: get_${NAME}_1
    procedure(get_v_i), deferred :: get_${NAME}_v
    generic, public              :: ${NAME} => get_${NAME}_1, get_${NAME}_v
  $endsub

  type, abstract :: therm_coeffs_t
   contains
     private
     $if($MPI)
     procedure(bcast_i), deferred, public :: bcast
     $endif
     $PROC_DECL(V_x2)
     $PROC_DECL(c_rad)
     $PROC_DECL(dc_rad)
     $PROC_DECL(c_gen)
     $PROC_DECL(c_thm)
     $PROC_DECL(nabla)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(dnabla_ad)
     $PROC_DECL(alpha_T)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(epsilon_ad)
     $PROC_DECL(epsilon_S)
  end type therm_coeffs_t

  ! Interfaces

  abstract interface

     subroutine bcast_i (this, root_rank)
       import therm_coeffs_t
       class(therm_coeffs_t), intent(inout) :: this
       integer, intent(in)                  :: root_rank
     end subroutine bcast_i

     function get_1_i (this, x) result (y)
       use core_kinds
       import therm_coeffs_t
       class(therm_coeffs_t), intent(in) :: this
       real(WP), intent(in)              :: x
       real(WP)                          :: y
     end function get_1_i

     function get_v_i (this, x) result (y)
       use core_kinds
       import therm_coeffs_t
       class(therm_coeffs_t), intent(in) :: this
       real(WP), intent(in)              :: x(:)
       real(WP)                          :: y(SIZE(x))
     end function get_v_i

  end interface

 ! Access specifiers

  private

  public :: therm_coeffs_t

end module gyre_therm_coeffs
