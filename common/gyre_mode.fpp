! Module   : gyre_mode
! Purpose  : mode classification
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

module gyre_mode

  ! Uses

  use core_kinds
  use core_parallel
  use core_hgroup

  use gyre_bvp
  use gyre_mech_coeffs
  use gyre_oscpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: mode_t
     private
     type(oscpar_t)           :: op
     real(WP), allocatable    :: x(:)
     complex(WP), allocatable :: xi_r(:)
     complex(WP), allocatable :: xi_h(:)
     complex(WP), allocatable :: phi_pri(:)
     complex(WP), allocatable :: dphi_pri(:)
     real(WP), allocatable    :: dE_dx(:)
     complex(WP), public      :: omega
     complex(WP), public      :: discrim
     real(WP), public         :: E
     integer                  :: n
     integer, public          :: n_p
     integer, public          :: n_g
   contains
     private
     procedure, public :: init
     procedure, public :: write
     procedure         :: write_gyre
  end type mode_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_md
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: mode_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init (this, mc, op, omega, discrim, x, xi_r, xi_h, phi_pri, dphi_pri)

    class(mode_t), intent(out)       :: this
    class(mech_coeffs_t), intent(in) :: mc
    type(oscpar_t), intent(in)       :: op
    complex(WP), intent(in)          :: omega
    complex(WP), intent(in)          :: discrim
    real(WP), intent(in)             :: x(:)
    complex(WP), intent(in)          :: xi_r(:)
    complex(WP), intent(in)          :: xi_h(:)
    complex(WP), intent(in)          :: phi_pri(:)
    complex(WP), intent(in)          :: dphi_pri(:)

    $CHECK_BOUNDS(SIZE(xi_r),SIZE(x))
    $CHECK_BOUNDS(SIZE(xi_h),SIZE(x))
    $CHECK_BOUNDS(SIZE(phi_pri),SIZE(x))
    $CHECK_BOUNDS(SIZE(dphi_pri),SIZE(x))

    ! Initialize the mode

    this%op = op

    this%xi_r = xi_r
    this%xi_h = xi_h
    this%phi_pri = phi_pri
    this%dphi_pri = dphi_pri

    this%omega = omega
    this%discrim = discrim

    this%dE_dx = kinetic(mc, this%op, omega, this%x, this%xi_r, this%xi_h)
    this%E = inertia(mc, this%op, omega, this%x, this%xi_r, this%xi_h, this%dE_dx)

    this%n = SIZE(this%x)

    call classify(this%x, REAL(this%xi_r), REAL(this%xi_h), this%n_p, this%n_g)

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_md (this, root_rank)

    class(mode_t), intent(inout) :: this
    integer, intent(in)          :: root_rank

    ! Broadcast the mode

    call bcast(this%op, root_rank)

    call bcast_alloc(this%x, root_rank)

    call bcast_alloc(this%xi_r, root_rank)
    call bcast_alloc(this%xi_h, root_rank)
    call bcast_alloc(this%phi_pri, root_rank)
    call bcast_alloc(this%dphi_pri, root_rank)

    call bcast_alloc(this%dE_dx, root_rank)

    call bcast(this%omega, root_rank)
    call bcast(this%discrim, root_rank)
    call bcast(this%E, root_rank)

    call bcast(this%n, root_rank)

    call bcast(this%n_p, root_rank)
    call bcast(this%n_g, root_rank)

  end subroutine bcast_md

  $endif

!****

  subroutine write (this, file, file_type)

    class(mode_t), intent(in)    :: this
    character(LEN=*), intent(in) :: file
    character(LEN=*), intent(in) :: file_type

    ! Write the mode

    select case (file_type)
    case ('GYRE')
       call write_gyre(this, file)
    case default
       $ABORT(Invalid file_type)
    end select

    ! Finish

    return

  end subroutine write

!****

  subroutine write_gyre (this, file)

    class(mode_t), intent(in)    :: this
    character(LEN=*), intent(in) :: file

    type(hgroup_t) :: hg

    ! Write the mode to a GYRE eigenfunction HDF5 file

    call hg%init(file, CREATE_FILE)

    call write_attr(hg, 'l', this%op%l)
    call write_attr(hg, 'lambda_0', this%op%lambda_0)

    call write_attr(hg, 'omega', this%omega)
    call write_attr(hg, 'discrim', this%discrim)

    call write_attr(hg, 'n', this%n)
    
    call write_attr(hg, 'n_p', this%n_p)
    call write_attr(hg, 'n_g', this%n_g)

    call write_dset(hg, 'x', this%x)

    call write_dset(hg, 'xi_r', this%xi_r)
    call write_dset(hg, 'xi_h', this%xi_h)
    call write_dset(hg, 'phi_pri', this%phi_pri)
    call write_dset(hg, 'dphi_pri', this%dphi_pri)

    call write_dset(hg, 'dE_dx', this%dE_dx)

    call hg%final()
    
    ! Finish

    return

  end subroutine write_gyre

!****

  subroutine classify (x, xi_r, xi_h, n_p, n_g)

    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: xi_r(:)
    real(WP), intent(in) :: xi_h(:)
    integer, intent(out) :: n_p
    integer, intent(out) :: n_g

    logical  :: inner_ext
    integer  :: i
    real(WP) :: y_2_cross

    $CHECK_BOUNDS(SIZE(xi_r),SIZE(x))
    $CHECK_BOUNDS(SIZE(xi_h),SIZE(x))

    ! Classify the non-radial eigenfunction using the
    ! Cowling-Scuflaire scheme

    n_p = 0
    n_g = 0
 
    inner_ext = ABS(xi_r(1)) > ABS(xi_r(2))

    x_loop : do i = 2,SIZE(x)-1

       ! If the innermost extremum in y_1 hasn't yet been reached,
       ! skip

       if(.NOT. inner_ext) then
          inner_ext = ABS(xi_r(i)) > ABS(xi_r(i-1)) .AND. ABS(xi_r(i)) > ABS(xi_r(i+1))
          cycle x_loop
       endif

       ! Look for a node in xi_r

       if(xi_r(i) >= 0._WP .AND. xi_r(i+1) < 0._WP) then

          y_2_cross = xi_h(i) - xi_r(i)*(xi_h(i+1) - xi_h(i))/(xi_r(i+1) - xi_r(i))

          if(y_2_cross >= 0._WP) then
             n_p = n_p + 1
          else
             n_g = n_g + 1
          endif

       elseif(xi_r(i) <= 0._WP .AND. xi_r(i+1) > 0._WP) then

         y_2_cross = xi_h(i) - xi_r(i)*(xi_h(i+1) - xi_h(i))/(xi_r(i+1) - xi_r(i))

          if(y_2_cross <= 0._WP) then
             n_p = n_p + 1
          else
             n_g = n_g + 1
          endif

       endif

    end do x_loop

    ! Finish

    return

  end subroutine classify

!*****

  function kinetic (mc, op, omega, x, xi_r, xi_h) result (dE_dx)

    class(mech_coeffs_t), intent(in) :: mc
    type(oscpar_t), intent(in)       :: op
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x(:)
    complex(WP), intent(in)          :: xi_r(:)
    complex(WP), intent(in)          :: xi_h(:)
    real(WP)                         :: dE_dx(SIZE(x))
    
    integer     :: i

    $CHECK_BOUNDS(SIZE(xi_r),SIZE(x))
    $CHECK_BOUNDS(SIZE(xi_h),SIZE(x))

    ! Calculate the kinetic energy density

    do i = 1,SIZE(x)

       associate(U => mc%U(x(i)), c_1 => mc%c_1(x(i)))

         dE_dx(i) = (ABS(xi_r(i))**2 + op%l*(op%l+1)*ABS(xi_h(i))**2)*U*x(i)**2/c_1

       end associate

    end do

    ! Finish

    return

  end function kinetic

!*****

  function inertia (mc, op, omega, x, xi_r, xi_h, dE_dx) result (E)

    class(mech_coeffs_t), intent(in) :: mc
    type(oscpar_t), intent(in)       :: op
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x(:)
    complex(WP), intent(in)          :: xi_r(:)
    complex(WP), intent(in)          :: xi_h(:)
    real(WP), intent(in)             :: dE_dx(:)
    real(WP)                         :: E

    integer     :: n
    real(WP)    :: E_norm

    $CHECK_BOUNDS(SIZE(xi_r),SIZE(x))
    $CHECK_BOUNDS(SIZE(xi_h),SIZE(x))

    $CHECK_BOUNDS(SIZE(dE_dx),SIZE(x))

    ! Calculate the normalized inertia using trapezoidal integration

    n = SIZE(x)

    E = SUM(0.5_WP*(dE_dx(2:) + dE_dx(:n-1))*(x(2:) - x(:n-1)))

    E_norm = ABS(xi_r(n))**2 + op%l*(op%l+1)*ABS(xi_h(n))**2

    if(E_norm == 0._WP) then
       $WARN(E_norm is zero, not normalizing inertia)
    else
       E = E/E_norm
    endif

    ! Finish

    return

  end function inertia

end module gyre_mode
