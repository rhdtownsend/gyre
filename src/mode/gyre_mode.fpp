! Module   : gyre_mode
! Purpose  : mode data
!
! Copyright 2013-2021 Rich Townsend & The GYRE Team
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
$include 'core_memory.inc'

module gyre_mode

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_context
  use gyre_constants
  use gyre_ext
  use gyre_grid
  use gyre_grid_util
  use gyre_math
  use gyre_mode_par
  use gyre_model
  use gyre_osc_par
  use gyre_point
  use gyre_state
  use gyre_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (wave_t) :: mode_t
     integer :: n_pg 
     integer :: n_p
     integer :: n_g
   contains
     private
     procedure :: classify_
  end type mode_t

  ! Interfaces

  interface mode_t
     module procedure mode_t_
  end interface mode_t

  interface reallocate
     module procedure reallocate_1_
  end interface reallocate

  ! Access specifiers

  private

  public :: mode_t
  public :: reallocate

  ! Procedures

contains

  function mode_t_ (wv) result (md)

    type(wave_t), intent(in) :: wv
    type(mode_t)             :: md

    complex(WP) :: y_1_ref
    complex(WP) :: f_phase

    ! Construct the mode_t

    md%wave_t = wv

    ! Normalize so that y_1 at the reference point is purely real, and
    ! the total inertia E is unity

    y_1_ref = md%y_i(1, md%j_ref)

    if (abs(y_1_ref) > TINY(0._WP)) then
       f_phase = CONJG(y_1_ref)/abs(y_1_ref)
    else
       f_phase = 1._WP
    endif

    md%scl = md%scl/sqrt(md%E())*f_phase

    ! Classify the mode

    call md%classify_()

    ! Finish

    return

  end function mode_t_

  !****

  subroutine classify_ (this)

    class(mode_t), intent(inout) :: this

    integer                 :: j
    real(WP)                :: y_1(this%n)
    real(WP)                :: y_2(this%n)
    integer                 :: j_i
    integer                 :: j_o
    real(WP)                :: x_i
    integer                 :: n_c
    integer                 :: n_a
    class(model_t), pointer :: ml
    type(grid_t)            :: gr

    ! Classify the mode based on its eigenfunctions

    if (this%l == 0) then

       ! Radial modes
       
       ! Look for the first monotonic segment in y_1 (this is to deal with
       ! noisy near-zero solutions at the origin)

       !$OMP PARALLEL DO
       do j = 1, this%n
          y_1(j) = REAL(this%y_i(1, j))
          y_2(j) = REAL(this%y_i(2, j))
       end do

       j_i = this%n

       mono_loop : do j = 2, this%n-1
          if ((y_1(j) >= y_1(j-1) .AND. y_1(j+1) >= y_1(j)) .OR. &
              (y_1(j) <= y_1(j-1) .AND. y_1(j+1) <= y_1(j))) then
             j_i = j
             exit mono_loop
          endif
       end do mono_loop

       ! Count winding numbers

       call count_windings_(y_1(j_i:), y_2(j_i:), n_c, n_a)

       ! Classify (the additional 1 is for the node at the center)

       this%n_p = n_a + n_c + 1
       this%n_g = 0

       this%n_pg = this%n_p - this%n_g

    elseif (this%l == 1 .AND. .NOT. this%os_p%alpha_grv == 0._WP) then

       ! Dipole modes (non-Cowling)

       ! Set up the Takata Y^a_1 and Y^a_2 functions

       !$OMP PARALLEL DO
       do j = 1, this%n
          y_1(j) = REAL(this%Yt_1(j))
          y_2(j) = REAL(this%Yt_2(j))
       end do

       ! Find the inner turning point (this is to deal with noisy
       ! near-zero solutions at the inner boundary)

       ml => this%model()

       call find_turn(this%context(), this%grid(), r_state_t(REAL(this%omega)), &
                      this%nm_p, this%os_p, j_i, x_i)

       ! Count winding numbers, taking care to avoid counting nodes at
       ! the center and surface

       if (y_1(this%n) == 0._WP) then
          j_o = this%n-1
       else
          j_o = this%n
       endif
       
       call count_windings_(y_1(j_i:j_o), y_2(j_i:j_o), n_c, n_a)

       ! Classify

       this%n_p = n_a
       this%n_g = n_c

       if (this%n_p >= this%n_g) then
          this%n_pg = this%n_p - this%n_g + 1
       else
          this%n_pg = this%n_p - this%n_g
       endif

    else

       ! Other modes

       !$OMP PARALLEL DO
       do j = 1, this%n
          y_1(j) = REAL(this%y_i(1, j))
          y_2(j) = REAL(this%y_i(2, j) + this%y_i(3, j))
       end do

       ! Handle special case where the inner boundary y_1 = 0 is
       ! appled off-center -- don't count the node there

       if (this%os_p%inner_bound == 'ZERO_R') then
          j_i = 2
       else
          j_i = 1
       endif

       j_o = this%n

       ! Count winding numbers

       call count_windings_(y_1(j_i:j_o), y_2(j_i:j_o), n_c, n_a)

       ! Classify

       this%n_p = n_a
       this%n_g = n_c

       this%n_pg = this%n_p - this%n_g

    endif

    ! Finish

    return

  contains

    subroutine count_windings_ (y_1, y_2, n_c, n_a, x)

      real(WP), intent(in)           :: y_1(:)
      real(WP), intent(in)           :: y_2(:)
      integer, intent(out)           :: n_c
      integer, intent(out)           :: n_a
      real(WP), optional, intent(in) :: x(:)

      integer  :: j
      real(WP) :: y_2_cross

      $CHECK_BOUNDS(SIZE(y_2),SIZE(y_1))

      if(PRESENT(x)) then
         $CHECK_BOUNDS(SIZE(x),SIZE(y_1))
      endif

      ! Count clockwise (n_c) and anticlockwise (n_a) windings in the (y_1,y_2) plane

      n_c = 0
      n_a = 0

      do j = 1,SIZE(y_1)-1

         ! Look for a node in y_1

         if (y_1(j) >= 0._WP .AND. y_1(j+1) < 0._WP) then

            ! Solve for the crossing ordinate

            y_2_cross = y_2(j) - y_1(j)*(y_2(j+1) - y_2(j))/(y_1(j+1) - y_1(j))

            if(y_2_cross >= 0._WP) then
               n_a = n_a + 1
               if(PRESENT(x)) print *,'A node:',x(j),x(j+1)
            else
               n_c = n_c + 1
               if(PRESENT(x)) print *,'C node:',x(j),x(j+1)
            endif

         elseif (y_1(j) <= 0._WP .AND. y_1(j+1) > 0._WP) then

            ! Solve for the crossing ordinate

            y_2_cross = y_2(j) - y_1(j)*(y_2(j+1) - y_2(j))/(y_1(j+1) - y_1(j))

            if (y_2_cross <= 0._WP) then
               n_a = n_a + 1
               if(PRESENT(x)) print *,'A node:',x(j),x(j+1)
            else
               n_c = n_c + 1
               if(PRESENT(x)) print *,'C node:',x(j),x(j+1)
            endif

         endif

      end do

      ! Finish

      return

    end subroutine count_windings_

  end subroutine classify_

  !****

  $REALLOCATE(type(mode_t),1)

end module gyre_mode
