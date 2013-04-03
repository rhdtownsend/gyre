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
     complex(WP), allocatable :: y(:,:)
     real(WP), allocatable    :: dE_dx(:)
     complex(WP), public      :: omega
     complex(WP), public      :: discrim
     real(WP), public         :: E
     integer                  :: n
     integer                  :: n_e
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

  subroutine init (this, mc, op, omega, discrim, x, y)

    class(mode_t), intent(out)       :: this
    class(mech_coeffs_t), intent(in) :: mc
    type(oscpar_t), intent(in)       :: op
    complex(WP), intent(in)          :: omega
    complex(WP), intent(in)          :: discrim
    real(WP), intent(in)             :: x(:)
    complex(WP), intent(in)          :: y(:,:)

    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Initialize the mode

    this%op = op

    this%x = x
    this%y = y

    this%omega = omega
    this%discrim = discrim

    this%dE_dx = kinetic(mc, this%op, omega, this%x, this%y(1:2,:))
    this%E = inertia(mc, this%op, omega, this%x, this%y(1:2,:), this%dE_dx)

    this%n = SIZE(this%x)
    this%n_e = SIZE(y, 1)

    call classify(this%x, REAL(this%y(1:2,:)), this%n_p, this%n_g)

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
    call bcast_alloc(this%y, root_rank)

    call bcast_alloc(this%dE_dx, root_rank)

    call bcast(this%omega, root_rank)
    call bcast(this%discrim, root_rank)
    call bcast(this%E, root_rank)

    call bcast(this%n, root_rank)
    call bcast(this%n_e, root_rank)

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
    call write_attr(hg, 'n_e', this%n_e)
    
    call write_attr(hg, 'n_p', this%n_p)
    call write_attr(hg, 'n_g', this%n_g)

    call write_dset(hg, 'x', this%x)
    call write_dset(hg, 'y', this%y)

    call write_dset(hg, 'dE_dx', this%dE_dx)

    call hg%final()
    
    ! Finish

    return

  end subroutine write_gyre

!****

  subroutine classify (x, y, n_p, n_g)

    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: y(:,:)
    integer, intent(out) :: n_p
    integer, intent(out) :: n_g

    logical  :: inner_ext
    integer  :: j
    real(WP) :: y_2_cross

    $CHECK_BOUNDS(SIZE(y, 1),2)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Classify the mode using the Cowling-Scuflaire scheme

    n_p = 0
    n_g = 0

    inner_ext = ABS(y(1,1)) > ABS(y(1,2))

    x_loop : do j = 2,SIZE(x)-1

       ! If the innermost extremum in y_1 hasn't yet been reached,
       ! skip

       if(.NOT. inner_ext) then
          inner_ext = ABS(y(1,j)) > ABS(y(1,j-1)) .AND. ABS(y(1,j)) > ABS(y(1,j+1))
          cycle x_loop
       endif

       ! Look for a node in y_1

       if(y(1,j) >= 0._WP .AND. y(1,j+1) < 0._WP) then

          y_2_cross = y(2,j) - y(1,j)*(y(2,j+1) - y(2,j))/(y(1,j+1) - y(1,j))

          if(y_2_cross >= 0._WP) then
             n_p = n_p + 1
          else
             n_g = n_g + 1
          endif

       elseif(y(1,j) <= 0._WP .AND. y(1,j+1) > 0._WP) then

         y_2_cross = y(2,j) - y(1,j)*(y(2,j+1) - y(2,j))/(y(1,j+1) - y(1,j))

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

  function kinetic (mc, op, omega, x, y) result (dE_dx)

    class(mech_coeffs_t), intent(in) :: mc
    type(oscpar_t), intent(in)       :: op
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x(:)
    complex(WP), intent(in)          :: y(:,:)
    real(WP)                         :: dE_dx(SIZE(x))
    
    integer     :: i
    real(WP)    :: U
    real(WP)    :: c_1
    complex(WP) :: xi_r
    complex(WP) :: xi_h

    $CHECK_BOUNDS(SIZE(y, 1),2)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Calculate the kinetic energy density

    !$OMP PARALLEL DO PRIVATE (U, c_1, xi_r, xi_h)
    do i = 1,SIZE(x)

       U = mc%U(x(i))
       c_1 = mc%c_1(x(i))

       if(op%l == 0) then
          xi_r = y(1,i)
          xi_h = 0._WP
       else
          xi_r = y(1,i)
          xi_h = y(2,i)/(c_1*omega**2)
       endif

       if(x(i) > 0._WP) then
          xi_r = xi_r*x(i)**(op%lambda_0+1._WP)
          xi_h = xi_h*x(i)**(op%lambda_0+1._WP)
       else
          if(op%lambda_0 /= -1._WP) then
             xi_r = 0._WP
             xi_h = 0._WP
          endif
       endif

       dE_dx(i) = (ABS(xi_r)**2 + op%l*(op%l+1)*ABS(xi_h)**2)*U*x(i)**2/c_1

    end do

    ! Finish

    return

  end function kinetic

!*****

  function inertia (mc, op, omega, x, y, dE_dx) result (E)

    class(mech_coeffs_t), intent(in) :: mc
    type(oscpar_t), intent(in)       :: op
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x(:)
    complex(WP), intent(in)          :: y(:,:)
    real(WP), intent(in)             :: dE_dx(:)
    real(WP)                         :: E

    integer     :: n
    complex(WP) :: xi_r
    complex(WP) :: xi_h
    real(WP)    :: E_norm

    $CHECK_BOUNDS(SIZE(y, 1),2)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    $CHECK_BOUNDS(SIZE(dE_dx),SIZE(x))

    ! Calculate the normalized inertia using trapezoidal integration

    n = SIZE(x)

    E = SUM(0.5_WP*(dE_dx(2:) + dE_dx(:n-1))*(x(2:) - x(:n-1)))

    if(op%l == 0) then
       xi_r = y(1,n)
       xi_h = 0._WP
    else
       xi_r = y(1,n)
       xi_h = y(2,n)/omega**2
    endif

    E_norm = ABS(xi_r)**2 + op%l*(op%l+1)*ABS(xi_h)**2

    if(E_norm == 0._WP) then
       $WARN(E_norm is zero, not normalizing inertia)
    else
       E = E/E_norm
    endif

    ! Finish

    return

  end function inertia

end module gyre_mode
