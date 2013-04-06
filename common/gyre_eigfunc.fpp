! Module   : gyre_eigfunc
! Purpose  : eigenfunction data
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

module gyre_eigfunc

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

  type :: eigfunc_t
     type(oscpar_t)           :: op
     real(WP), allocatable    :: x(:)
     complex(WP), allocatable :: xi_r(:)
     complex(WP), allocatable :: xi_h(:)
     complex(WP), allocatable :: phi_pri(:)
     complex(WP), allocatable :: dphi_pri(:)
     complex(WP), allocatable :: del_S(:)
     complex(WP), allocatable :: del_L(:)
     complex(WP)              :: omega
     integer                  :: n
   contains
     private
     procedure, public :: init
     procedure, public :: write
     procedure         :: write_gyre
     procedure, public :: classify
     procedure, public :: kinetic
     procedure, public :: inertia
  end type eigfunc_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_ef
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: eigfunc_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init (this, op, omega, x, xi_r, xi_h, phi_pri, dphi_pri, del_S, del_L)

    class(eigfunc_t), intent(out)    :: this
    type(oscpar_t), intent(in)       :: op
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x(:)
    complex(WP), intent(in)          :: xi_r(:)
    complex(WP), intent(in)          :: xi_h(:)
    complex(WP), intent(in)          :: phi_pri(:)
    complex(WP), intent(in)          :: dphi_pri(:)
    complex(WP), intent(in)          :: del_S(:)
    complex(WP), intent(in)          :: del_L(:)

    $CHECK_BOUNDS(SIZE(xi_r),SIZE(x))
    $CHECK_BOUNDS(SIZE(xi_h),SIZE(x))
    $CHECK_BOUNDS(SIZE(phi_pri),SIZE(x))
    $CHECK_BOUNDS(SIZE(dphi_pri),SIZE(x))
    $CHECK_BOUNDS(SIZE(del_S),SIZE(x))
    $CHECK_BOUNDS(SIZE(del_L),SIZE(x))

    ! Initialize the eigfunc

    this%op = op

    this%x = x

    this%xi_r = xi_r
    this%xi_h = xi_h
    this%phi_pri = phi_pri
    this%dphi_pri = dphi_pri
    this%del_S = del_S
    this%del_L = del_L

    this%omega = omega

    this%n = SIZE(this%x)

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_ef (this, root_rank)

    class(eigfunc_t), intent(inout) :: this
    integer, intent(in)             :: root_rank

    ! Broadcast the eigfunc

    call bcast(this%op, root_rank)

    call bcast_alloc(this%x, root_rank)

    call bcast_alloc(this%xi_r, root_rank)
    call bcast_alloc(this%xi_h, root_rank)
    call bcast_alloc(this%phi_pri, root_rank)
    call bcast_alloc(this%dphi_pri, root_rank)
    call bcast_alloc(this%del_S, root_rank)
    call bcast_alloc(this%del_L, root_rank)

    call bcast(this%n, root_rank)

  end subroutine bcast_ef

  $endif

!****

  subroutine write (this, file, file_type)

    class(eigfunc_t), intent(in) :: this
    character(LEN=*), intent(in) :: file
    character(LEN=*), intent(in) :: file_type

    ! Write the eigfunc

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

    class(eigfunc_t), intent(in) :: this
    character(LEN=*), intent(in) :: file

    type(hgroup_t) :: hg

    ! Write the eigfunc to a GYRE eigenfunction HDF5 file

    call hg%init(file, CREATE_FILE)

    call write_attr(hg, 'l', this%op%l)

    call write_attr(hg, 'omega', this%omega)

    call write_attr(hg, 'n', this%n)
    
    call write_dset(hg, 'x', this%x)

    call write_dset(hg, 'xi_r', this%xi_r)
    call write_dset(hg, 'xi_h', this%xi_h)
    call write_dset(hg, 'phi_pri', this%phi_pri)
    call write_dset(hg, 'dphi_pri', this%dphi_pri)
    call write_dset(hg, 'del_S', this%del_S)
    call write_dset(hg, 'del_L', this%del_L)

    call hg%final()
    
    ! Finish

    return

  end subroutine write_gyre

!****

  subroutine classify (this, n_p, n_g)

    class(eigfunc_t), intent(in) :: this
    integer, intent(out)         :: n_p
    integer, intent(out)         :: n_g

    real(WP) :: xi_r(this%n)
    real(WP) :: xi_h(this%n)
    logical  :: inner_ext
    integer  :: i
    real(WP) :: y_2_cross

    ! Classify the eigenfunction using the Cowling-Scuflaire scheme

    xi_r = REAL(this%xi_r)
    xi_h = REAL(this%xi_h)

    n_p = 0
    n_g = 0
 
    inner_ext = ABS(xi_r(1)) > ABS(this%xi_r(2))

    x_loop : do i = 2,this%n-1

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

  function kinetic (this, mc) result (dE_dx)

    class(eigfunc_t), intent(in)     :: this
    class(mech_coeffs_t), intent(in) :: mc
    real(WP)                         :: dE_dx(this%n)
    
    integer     :: i

    ! Calculate the kinetic energy density

    do i = 1,this%n
       associate(U => mc%U(this%x(i)), c_1 => mc%c_1(this%x(i)))
         dE_dx(i) = (ABS(this%xi_r(i))**2 + this%op%l*(this%op%l+1)*ABS(this%xi_h(i))**2)*U*this%x(i)**2/c_1
       end associate
    end do

    ! Finish

    return

  end function kinetic

!*****

  function inertia (this, mc) result (E)

    class(eigfunc_t), intent(in)     :: this
    class(mech_coeffs_t), intent(in) :: mc
    real(WP)                         :: E

    real(WP) :: dE_dx(this%n)
    real(WP) :: E_norm

    ! Calculate the kinetic energy density

    dE_dx = this%kinetic(mc)

    ! Integrate it to obtain the mode inertia

    E = SUM(0.5_WP*(dE_dx(2:) + dE_dx(:this%n-1))*(this%x(2:) - this%x(:this%n-1)))

    ! Normalize

    E_norm = ABS(this%xi_r(this%n))**2 + this%op%l*(this%op%l+1)*ABS(this%xi_h(this%n))**2

    if(E_norm == 0._WP) then
       $WARN(E_norm is zero, not normalizing inertia)
    else
       E = E/E_norm
    endif

    ! Finish

    return

  end function inertia

end module gyre_eigfunc
