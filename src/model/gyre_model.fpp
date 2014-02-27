! Module   : gyre_model
! Purpose  : stellar model (interface)
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

module gyre_model

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure(y_1_), deferred :: ${NAME}_1_
    procedure(y_v_), deferred :: ${NAME}_v_
    generic, public           :: ${NAME} => ${NAME}_1_, ${NAME}_v_
  $endsub

  type, abstract :: model_t
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
     procedure(pi_c_), deferred, public         :: pi_c
     procedure(is_zero_), deferred, public      :: is_zero
     procedure(attach_cache_), deferred, public :: attach_cache
     procedure(detach_cache_), deferred, public :: detach_cache
     procedure(fill_cache_), deferred, public   :: fill_cache
     procedure, public                          :: omega => omega_
     procedure                                  :: omega_c_1_
     procedure                                  :: omega_c_v_
     generic, public                            :: omega_c => omega_c_1_, omega_c_v_
  end type model_t

  ! Interfaces

  abstract interface

     function y_1_ (this, x) result (y)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       real(WP), intent(in)       :: x
       real(WP)                   :: y
     end function y_1_

     function y_v_ (this, x) result (y)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       real(WP), intent(in)       :: x(:)
       real(WP)                   :: y(SIZE(x))
     end function y_v_

     function pi_c_ (this) result (pi_c)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       real(WP)                   :: pi_c
     end function pi_c_

     function is_zero_ (this, x) result (is_zero)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       real(WP), intent(in)       :: x
       logical                    :: is_zero
     end function is_zero_

     subroutine attach_cache_ (this, cc)
       use gyre_cocache
       import model_t
       class(model_t), intent(inout)         :: this
       class(cocache_t), pointer, intent(in) :: cc
     end subroutine attach_cache_

     subroutine detach_cache_ (this)
       import model_t
       class(model_t), intent(inout) :: this
     end subroutine detach_cache_

     subroutine fill_cache_ (this, x)
       use core_kinds
       import model_t
       class(model_t), intent(inout) :: this
       real(WP), intent(in)          :: x(:)
     end subroutine fill_cache_

  end interface

 ! Access specifiers

  private

  public :: model_t

  ! Procedures

contains

  $if($GFORTRAN_PR57922)

  subroutine final (this)

    class(model_t), intent(inout) :: this

    ! Finalize the model_t

    ! Finish

    return

  end subroutine final

  $endif

!****

  function omega_ (this, x, m, omega_c) result (omega)

    class(model_t), intent(in) :: this
    real(WP), intent(in)       :: x
    integer, intent(in)        :: m
    complex(WP), intent(in)    :: omega_c
    complex(WP)                :: omega

    ! Calculate the intertial frequency from the co-rotating frequency

    omega = omega_c + m*this%Omega_rot(x)

    ! Finish

    return

  end function omega_

!****

  function omega_c_1_ (this, x, m, omega) result (omega_c)

    class(model_t), intent(in) :: this
    real(WP), intent(in)       :: x
    integer, intent(in)        :: m
    complex(WP), intent(in)    :: omega
    complex(WP)                :: omega_c

    ! Calculate the co-rotating frequency from the inertial frequency

    omega_c = omega - m*this%Omega_rot(x)

    ! Finish

    return

  end function omega_c_1_

!****

  function omega_c_v_ (this, x, m, omega) result (omega_c)

    class(model_t), intent(in) :: this
    real(WP), intent(in)       :: x(:)
    integer, intent(in)        :: m
    complex(WP), intent(in)    :: omega
    complex(WP)                :: omega_c(SIZE(x))

    ! Calculate the co-rotating frequency from the inertial frequency

    omega_c = omega - m*this%Omega_rot(x)

    ! Finish

    return

  end function omega_c_v_

end module gyre_model
