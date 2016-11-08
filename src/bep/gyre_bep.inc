! Incfile  : gyre_bep
! Purpose  : boundary eigenvalue problems (template)
!
! Copyright 2013-2015 Rich Townsend
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_bep

  ! Uses

  use core_kinds

  use gyre_bound
  use gyre_diff
  use gyre_ext
  use gyre_num_par
  use gyre_status
  use gyre_sysmtx
  use gyre_sysmtx_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ${T}_bep_t
     private
     class(${T}_diff_t), allocatable   :: df(:)
     class(${T}_bound_t), allocatable  :: bd
     class(${T}_sysmtx_t), allocatable :: sm
     real(WP)                          :: omega_min
     real(WP)                          :: omega_max
     integer, public                   :: n_k
     integer, public                   :: n_e
     integer, public                   :: n_i
     integer, public                   :: n_o
   contains 
     private
     procedure, public :: eval_discrim
     procedure         :: build_
     procedure, public :: solve
     procedure, public :: solve_forced
  end type ${T}_bep_t

  ! Interfaces

  interface ${T}_bep_t
     module procedure ${T}_bep_t
  end interface ${T}_bep_t

  ! Access specifiers

  private

  public :: ${T}_bep_t
  public :: ${T}_bep_t_

contains

  function ${T}_bep_t_ (bd, df, omega_min, omega_max, nm_p) result (bp)

    class(${T}_diff_t), intent(in)   :: df(:)
    class(${T}_bound_t), intent(in)  :: bd
    real(WP), intent(in)             :: omega_min
    real(WP), intent(in)             :: omega_max
    type(num_par_t), intent(in)      :: nm_p
    type(${T}_bep_t)                 :: bp

    integer :: n_k
    integer :: n_e
    integer :: n_i
    integer :: n_o

    ! Perform basic validations

    n_k = SIZE(df) + 1
    n_e = bd%n_e

    n_i = bd%n_i
    n_o = bd%n_o

    $CHECK_BOUNDS(n_i+n_o,n_e)

    $ASSERT_DEBUG(ALL(df%n_e == n_e),Dimension mismatch)

    ! Construct the bep_t

    allocate(bp%df(n_k-1), SOURCE=df)

    allocate(bp%bd, SOURCE=bd)

    allocate(bp%sm, SOURCE=${T}_sysmtx_t(n_k-1, n_e, n_i, n_o, nm_p))

    bp%n_k = n_k
    bp%n_e = n_e
    bp%n_i = n_i
    bp%n_o = n_o

    bp%omega_min = omega_min
    bp%omega_max = omega_max
    
    ! Finish
    
    return

  end function ${T}_bep_t_

  !****

  subroutine eval_discrim (this, omega, discrim, status)

    class(${T}_bep_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: omega
    type(${T}_ext_t), intent(out)    :: discrim
    integer, intent(out)             :: status

    ! Evaluate the discriminant as the determinant of the sysmtx

    if (REAL(omega, WP) >= this%omega_min .AND. REAL(omega, WP) <= this%omega_max) then

       call this%build_(omega)

       call this%sm%factorize()
       discrim = this%sm%det()

       status = STATUS_OK

    else

       status = STATUS_OMEGA_DOMAIN

    endif

    ! Finish

    return

  end subroutine eval_discrim

  !****

  subroutine build_ (this, omega)

    class(${T}_bep_t), target, intent(inout) :: this
    $TYPE(WP), intent(in)                    :: omega

    $TYPE(WP)        :: B_i(this%n_i,this%n_e)
    $TYPE(WP)        :: B_o(this%n_o,this%n_e)
    $TYPE(WP)        :: E_l(this%n_e,this%n_e)
    $TYPE(WP)        :: E_r(this%n_e,this%n_e)
    type(${T}_ext_t) :: scl
    integer          :: k
    
    ! Build the sysmtx

    ! Boundary conditions

    call this%bd%build_i(omega, B_i, scl)
    call this%sm%set_B_i(B_i, scl)

    call this%bd%build_o(omega, B_o, scl)
    call this%sm%set_B_o(B_o, scl)

    ! Difference equations

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, scl) SCHEDULE (DYNAMIC)
    sub_loop : do k = 1, this%n_k-1
       call this%df(k)%build(omega, E_l, E_r, scl)
       call this%sm%set_E(k, E_l, E_r, scl)
    end do sub_loop

    ! Finish

    return

  end subroutine build_

  !****

  subroutine solve (this, omega, y, discrim)

    class(${T}_bep_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: omega
    $TYPE(WP), intent(out)           :: y(:,:)
    type(${T}_ext_t), intent(out)    :: discrim

    $TYPE(WP) :: v(this%n_e*this%n_k)

    $CHECK_BOUNDS(SIZE(y, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),this%n_k)

    ! Solve for the null vector of the sysmtx

    call this%build_(omega)

    call this%sm%factorize()
    discrim = this%sm%det()
    v = this%sm%null_vec()

    ! Pack the null vector into the solution vector

    y = RESHAPE(v, SHAPE(y))

    ! Finish

    return

  end subroutine solve

  !****

  subroutine solve_forced (this, omega, w_nz, y)

    class(${T}_bep_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: omega
    $TYPE(WP), intent(in)            :: w_nz(:)
    $TYPE(WP), intent(out)           :: y(:,:)

    $TYPE(WP) :: v(this%n_e*this%n_k)

    $CHECK_BOUNDS(SIZE(w_nz),this%n_e)

    $CHECK_BOUNDS(SIZE(y, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),this%n_k)

    ! Solve for the solution vector

    call this%build_(omega)

    call this%sm%factorize()

    v = this%sm%soln_vec(w_nz)

    ! Pack the solution vector into the solution vector

    y = RESHAPE(v, SHAPE(y))
    
    ! Finish

    return

  end subroutine solve_forced

end module gyre_${T}_bep

$endsub
