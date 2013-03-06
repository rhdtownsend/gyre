! Module   : gyre_bvp
! Purpose  : solve boundary-value problems (interface)

$include 'core.inc'

module gyre_bvp

  ! Uses

  use core_kinds

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: bvp_t
   contains 
     private
     procedure(discrim_i), deferred, public :: discrim
     procedure(recon_i), deferred, public   :: recon
  end type bvp_t

  ! Interfaces

  abstract interface

     subroutine recon_i (this, omega, x, y)
       use core_kinds
       import bvp_t
       class(bvp_t), intent(inout)           :: this
       complex(WP), intent(in)               :: omega
       real(WP), allocatable, intent(out)    :: x(:)
       complex(WP), allocatable, intent(out) :: y(:,:)
     end subroutine recon_i
       
     function discrim_i (this, omega) result (discrim)
       use core_kinds
       use gyre_ext_arith
       import bvp_t
       class(bvp_t), intent(inout) :: this
       complex(WP), intent(in)     :: omega
       type(ext_complex_t)         :: discrim
     end function discrim_i

  end interface

  ! Access specifiers

  private

  public :: bvp_t

end module gyre_bvp
