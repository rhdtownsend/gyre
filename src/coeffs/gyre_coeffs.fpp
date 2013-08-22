! Module   : gyre_coeffs
! Purpose  : structure coefficients (interface)
!
! Copyright 2013 Rich Townsend
!
! This file is part of GYRE. GYRE is free software: you can
! redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, version 3.
!
! GYRE is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

$include 'core.inc'

module gyre_coeffs

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure(y_1_i), deferred :: ${NAME}_1
    procedure(y_v_i), deferred :: ${NAME}_v
    generic, public            :: ${NAME} => ${NAME}_1, ${NAME}_v
  $endsub

  type, abstract :: coeffs_t
   contains
     private
     $if($GFORTRAN_PR57922)
     procedure, public :: final
     $endif
     $PROC_DECL(V)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(delta)
     $PROC_DECL(c_rad)
     $PROC_DECL(dc_rad)
     $PROC_DECL(c_thm)
     $PROC_DECL(c_dif)
     $PROC_DECL(c_eps_ad)
     $PROC_DECL(c_eps_S)
     $PROC_DECL(nabla)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(tau_thm)
     $PROC_DECL(Omega_rot)
     procedure(pi_c_i), deferred, public         :: pi_c
     procedure(enable_cache_i), deferred, public :: enable_cache
     procedure(enable_cache_i), deferred, public :: disable_cache
     procedure(fill_cache_i), deferred, public   :: fill_cache
     procedure, public                           :: omega
     procedure, public                           :: omega_c
  end type coeffs_t

  ! Interfaces

  abstract interface

     function y_1_i (this, x) result (y)
       use core_kinds
       import coeffs_t
       class(coeffs_t), intent(in) :: this
       real(WP), intent(in)        :: x
       real(WP)                    :: y
     end function y_1_i

     function y_v_i (this, x) result (y)
       use core_kinds
       import coeffs_t
       class(coeffs_t), intent(in) :: this
       real(WP), intent(in)        :: x(:)
       real(WP)                    :: y(SIZE(x))
     end function y_v_i

     function pi_c_i (this) result (pi_c)
       use core_kinds
       import coeffs_t
       class(coeffs_t), intent(in) :: this
       real(WP)                    :: pi_c
     end function pi_c_i

     subroutine enable_cache_i (this)
       import coeffs_t
       class(coeffs_t), intent(inout) :: this
     end subroutine enable_cache_i

     subroutine fill_cache_i (this, x)
       use core_kinds
       import coeffs_t
       class(coeffs_t), intent(inout) :: this
       real(WP), intent(in)           :: x(:)
     end subroutine fill_cache_i

  end interface

 ! Access specifiers

  private

  public :: coeffs_t

  ! Procedures

contains

  $if($GFORTRAN_PR57922)

  subroutine final (this)

    class(coeffs_t), intent(inout) :: this

    ! Finalize the coeffs

    ! Finish

    return

  end subroutine final

  $endif

!****

  function omega (this, x, m, omega_c)

    class(coeffs_t), intent(in) :: this
    real(WP), intent(in)        :: x
    integer, intent(in)         :: m
    complex(WP), intent(in)     :: omega_c
    complex(WP)                 :: omega

    ! Calculate the intertial frequency from the co-rotating frequency

    omega = omega_c + m*this%Omega_rot(x)

    ! Finish

    return

  end function omega

!****

  function omega_c (this, x, m, omega)

    class(coeffs_t), intent(in) :: this
    real(WP), intent(in)        :: x
    integer, intent(in)         :: m
    complex(WP), intent(in)     :: omega
    complex(WP)                 :: omega_c

    ! Calculate the co-rotating frequency from the inertial frequency

    omega_c = omega - m*this%Omega_rot(x)

    ! Finish

    return

  end function omega_c

end module gyre_coeffs
