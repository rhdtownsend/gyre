! Module   : gyre_m_discrim_func
! Purpose  : modulus of complex discriminant function
!
! Copyright 2018-2020 Rich Townsend & The GYRE Team
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

module gyre_m_discrim_func

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_ext
  use gyre_ext_func
  use gyre_math
  use gyre_status
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_ext_func_t) :: m_discrim_func_t
     private
     class(c_bvp_t), pointer       :: bp
     class(c_state_t), allocatable :: st
     real(WP)                      :: omega_min
     real(WP)                      :: omega_max
   contains 
     private
     procedure, public :: eval => eval_
  end type m_discrim_func_t

  ! Interfaces

  interface m_discrim_func_t
     module procedure m_discrim_func_t_
  end interface m_discrim_func_t

  ! Access specifiers

  private

  public :: m_discrim_func_t

  ! Procedures

contains

  function m_discrim_func_t_ (bp, st, omega_min, omega_max) result (df)

    class(c_bvp_t), pointer, intent(in) :: bp
    class(c_state_t), intent(in)        :: st
    real(WP), intent(in)                :: omega_min
    real(WP), intent(in)                :: omega_max
    type(m_discrim_func_t)              :: df

    ! Construct the m_discrim_func_t

    df%bp => bp

    allocate(df%st, SOURCE=st)

    df%omega_min = omega_min
    df%omega_max = omega_max

    ! Finish

    return

  end function m_discrim_func_t_

  !****

  subroutine eval_ (this, rx, f_rx, status)

    class(m_discrim_func_t), intent(inout) :: this
    type(r_ext_t), intent(in)              :: rx
    type(r_ext_t), intent(out)             :: f_rx
    integer, intent(out)                   :: status

    complex(WP)     :: omega
    type(c_state_t) :: st

    ! Evaluate the discriminant function

    omega = CMPLX(real(rx), KIND=WP)

    if (REAL(omega) >= this%omega_min .AND. REAL(omega) <= this%omega_max) then

       this%st%omega = omega
       this%st%omega_r = REAL(omega)

       call this%bp%build(this%st)
       call this%bp%factor()

       f_rx = abs(this%bp%det())

       status = STATUS_OK

    else

       status = STATUS_OMEGA_DOMAIN

    endif
    
    ! Finish

    return

  end subroutine eval_

end module gyre_m_discrim_func
