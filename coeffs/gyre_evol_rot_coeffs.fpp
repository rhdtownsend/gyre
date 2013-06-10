! Module   : gyre_evol_rot_coeffs
! Purpose  : rotation coefficients for evolutionary models
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

module gyre_evol_rot_coeffs

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  use core_spline

  use gyre_rot_coeffs
  use gyre_cocache

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $VAR_DECL $sub
    $local $NAME $1
    type(spline_t) :: sp_$NAME
  $endsub
  
  $define $PROC_DECL $sub
    $local $NAME $1
    procedure :: ${NAME}_1
    procedure :: ${NAME}_v
  $endsub

  $define $PROC_DECL_GEN $sub
    $local $NAME $1
    procedure       :: ${NAME}_1
    procedure       :: ${NAME}_v
    generic, public :: ${NAME} => ${NAME}_1, ${NAME}_v
  $endsub

  type, extends(rot_coeffs_t) :: evol_rot_coeffs_t
     private
     type(cocache_t) :: cc
     $VAR_DECL(Omega_rot)
     real(WP), public :: M_star
     real(WP), public :: R_star
     real(WP), public :: G
     logical          :: cc_enabled
   contains
     private
     procedure, public :: init
     $PROC_DECL(Omega_rot)
     procedure, public :: enable_cache
     procedure, public :: disable_cache
     procedure, public :: fill_cache
  end type evol_rot_coeffs_t
 
  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_rc
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: evol_rot_coeffs_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  recursive subroutine init (this, G, M_star, R_star, r, Omega_rot, deriv_type, add_center)

    class(evol_rot_coeffs_t), intent(out) :: this
    real(WP), intent(in)                  :: G
    real(WP), intent(in)                  :: M_star
    real(WP), intent(in)                  :: R_star
    real(WP), intent(in)                  :: r(:)
    real(WP), intent(in)                  :: Omega_rot(:)
    character(LEN=*), intent(in)          :: deriv_type
    logical, intent(in), optional         :: add_center

    logical  :: add_center_
    real(WP) :: Omega_rot_(SIZE(r))
    real(WP) :: x(SIZE(r))

    $CHECK_BOUNDS(SIZE(Omega_rot),SIZE(r))

    if(PRESENT(add_center)) then
       add_center_ = add_center
    else
       add_center_ = .FALSE.
    endif

    ! See if we need a central point

    if(add_center_) then

       ! Add a central point and initialize using recursion

       call this%init(G, M_star, R_star, [0._WP,r], y_centered(r, Omega_rot), deriv_type, .FALSE.)

    else

       ! Perform basic validations
       
       $ASSERT(r(1) == 0._WP,First grid point not at center)

       $ASSERT(ALL(r(2:) >= r(:SIZE(r)-1)),Non-monotonic radius data)

       ! Calculate the dimensionless rotation frequency

       Omega_rot_ = SQRT(R_star**3/(G*M_star))*Omega_rot

       x = r/R_star

       ! Initialize the rot_coeffs

       !$OMP PARALLEL SECTIONS
       !$OMP SECTION
       call this%sp_Omega_rot%init(x, Omega_rot_, deriv_type, dy_dx_a=0._WP)
       !$OMP END PARALLEL SECTIONS

       this%M_star = M_star
       this%R_star = R_star

       this%G = G

       this%cc_enabled = .FALSE.

    endif

    ! Finish

    return

  contains

    function y_centered (x, y)
      
      real(WP), intent(in) :: x(:)
      real(WP), intent(in) :: y(:)
      real(WP)             :: y_centered(SIZE(y)+1)

      real(WP) :: y_center

      $CHECK_BOUNDS(SIZE(x),SIZE(y))

      $ASSERT(SIZE(y) >= 2,Insufficient grid points)

      ! Use parabola fitting to interpolate y at the center
      
      y_center = (x(2)**2*y(1) - x(1)**2*y(2))/(x(2)**2 - x(1)**2)

      ! Create the centered array

      y_centered = [y_center,y]

      ! Finish

      return

    end function y_centered

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_rc (rc, root_rank)

    class(evol_rot_coeffs_t), intent(inout) :: rc
    integer, intent(in)                     :: root_rank

    ! Broadcast the rot_coeffs

    call bcast(rc%cc, root_rank)

    call bcast(rc%sp_Omega_rot, root_rank)

    call bcast(rc%M_star, root_rank)
    call bcast(rc%R_star, root_rank)

    call bcast(rc%G, root_rank)

    call bcast(rc%cc_enabled, root_rank)

    ! Finish

    return

  end subroutine bcast_rc

  $endif

!****

  $define $PROC $sub

  $local $NAME $1
  $local $I_CC $2

  function ${NAME}_1 (this, x) result ($NAME)

    class(evol_rot_coeffs_t), intent(in) :: this
    real(WP), intent(in)                 :: x
    real(WP)                             :: $NAME

    ! Interpolate $NAME

    if(this%cc_enabled) then
       $NAME = this%cc%lookup($I_CC, x)
    else
       $NAME = this%sp_$NAME%interp(x)
    endif
       
    ! Finish

    return

  end function ${NAME}_1

!****

  function ${NAME}_v (this, x) result ($NAME)

    class(evol_rot_coeffs_t), intent(in) :: this
    real(WP), intent(in)                 :: x(:)
    real(WP)                             :: $NAME(SIZE(x))

    ! Interpolate $NAME

    $NAME = this%sp_$NAME%interp(x)

    ! Finish

    return

  end function ${NAME}_v

  $endsub

  $PROC(Omega_rot,1)

!****

  subroutine enable_cache (this)

    class(evol_rot_coeffs_t), intent(inout) :: this

    ! Enable the coefficient cache

    this%cc_enabled = .TRUE.

    ! Finish

    return

  end subroutine enable_cache

!****

  subroutine disable_cache (this)

    class(evol_rot_coeffs_t), intent(inout) :: this

    ! Disable the coefficient cache

    this%cc_enabled = .FALSE.

    ! Finish

    return

  end subroutine disable_cache

!****

  subroutine fill_cache (this, x)

    class(evol_rot_coeffs_t), intent(inout) :: this
    real(WP), intent(in)                    :: x(:)

    real(WP) :: c(11,SIZE(x))

    ! Fill the coefficient cache

    !$OMP PARALLEL SECTIONS
    !$OMP SECTION
    c(1,:) = this%Omega_rot(x)
    !$OMP END PARALLEL SECTIONS

    call this%cc%init(x, c)

    ! Finish

    return

  end subroutine fill_cache

end module gyre_evol_rot_coeffs
