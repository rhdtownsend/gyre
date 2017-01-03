! Module   : gyre_ad_bvp
! Purpose  : adiabatic bounary value problem solver
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_ad_bvp

  ! Uses

  use core_kinds

  use gyre_ad_bound
  use gyre_ad_eqns
  use gyre_ad_diff
  use gyre_ad_vars
  use gyre_bvp
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_model
  use gyre_mode
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_soln

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_bvp_t) :: ad_bvp_t
     class(model_t), pointer :: ml => null()
     type(grid_t)            :: gr
     type(ad_eqns_t)         :: eq
     type(ad_vars_t)         :: vr
     type(mode_par_t)        :: md_p
     type(osc_par_t)         :: os_p
  end type ad_bvp_t

  ! Interfaces

  interface ad_bvp_t
     module procedure ad_bvp_t_
  end interface ad_bvp_t

  interface mode_t
     module procedure mode_t_
  end interface mode_t

  ! Access specifiers

  private

  public :: ad_bvp_t
  public :: mode_t

  ! Procedures

contains

  function ad_bvp_t_ (ml, gr, md_p, nm_p, os_p) result (bp)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    type(mode_par_t), intent(in)        :: md_p
    type(num_par_t), intent(in)         :: nm_p
    type(osc_par_t), intent(in)         :: os_p
    type(ad_bvp_t)                      :: bp

    type(ad_bound_t)             :: bd
    integer                      :: k
    type(ad_diff_t), allocatable :: df(:)
    real(WP)                     :: omega_min
    real(WP)                     :: omega_max

    ! Construct the ad_bvp_t

    ! Initialize the boundary conditions

    bd = ad_bound_t(ml, gr, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_k-1))

    do k = 1, gr%n_k-1
       df(k) = ad_diff_t(ml, gr, k, md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t_(bd, df, nm_p) 

    ! Other initializations
 
    bp%ml => ml
    bp%gr = gr

    bp%eq = ad_eqns_t(ml, gr, md_p, os_p)
    bp%vr = ad_vars_t(ml, gr, md_p, os_p)

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function ad_bvp_t_

  !****

    ! if (REAL(omega, WP) >= this%omega_min .AND. REAL(omega, WP) <= this%omega_max) then

    !    call this%build_(omega)

    !    call this%sm%factorize()
    !    discrim = this%sm%det()

    !    status = STATUS_OK

    ! else

    !    status = STATUS_OMEGA_DOMAIN

    ! endif

  !****

  function mode_t_ (bp, j, omega) result (md)

    class(ad_bvp_t), intent(inout) :: bp
    integer, intent(in)            :: j
    real(WP), intent(in)           :: omega
    type(mode_t)                   :: md

    type(soln_t) :: sl

    ! Construct the mode_t

    sl = soln_t_(bp, omega)

    md = mode_t(bp%ml, sl, j, bp%md_p, bp%os_p)

    ! Finish

    return

  end function mode_t_

  !****

  function soln_t_ (bp, omega) result (sl)

    class(ad_bvp_t), intent(inout) :: bp
    real(WP), intent(in)           :: omega
    type(soln_t)                   :: sl

    real(WP)      :: y(4,bp%n_k)
    type(r_ext_t) :: discrim
    integer       :: k
    real(WP)      :: xA(4,4)
    real(WP)      :: dy_dx(4,bp%n_k)
    real(WP)      :: H(4,4)
    real(WP)      :: dH(4,4)
    complex(WP)   :: y_c(4,bp%n_k)
    complex(WP)   :: dy_c_dx(4,bp%n_k)
    integer       :: i

    ! Calculate the solution vector and discriminant

    call bp%build(omega)

    y = bp%soln_vec()
    discrim = bp%det()

    ! Calculate its derivatives (nb: the vacuum check prevents errors, but
    ! leads to incorrect values for dy_dx)

    !$OMP PARALLEL DO PRIVATE (xA)
    do k = 1, bp%n_k

       associate (pt => bp%gr%pt(k))

         if (bp%ml%vacuum(pt)) then
            xA = 0._WP
         else
            xA = bp%eq%xA(pt, omega)
         endif

         if (pt%x /= 0._WP) then
            dy_dx(:,k) = MATMUL(xA, y(:,k))/pt%x
         else
            dy_dx(:,k) = 0._WP
         endif

       end associate

    end do

    ! Convert to canonical form (nb: the vacuum check prevents errors, but
    ! leads to incorrect values for dy_c_dx)

    !$OMP PARALLEL DO PRIVATE (H, dH)
    do k = 1, bp%n_k

       associate (pt => bp%gr%pt(k))

         H = bp%vr%H(pt, omega)

         y_c(:,k) = MATMUL(H, y(:,k))

         if (bp%ml%vacuum(pt)) then
            dH = 0._WP
         else
            dH = bp%vr%dH(pt, omega)
         endif

         if (pt%x /= 0._WP) then
            dy_c_dx(:,k) = MATMUL(dH/pt%x, y(:,k)) + MATMUL(H, dy_dx(:,k))
         else
            dy_c_dx(:,k) = 0._WP
         endif

       end associate

    end do

    ! Construct the soln_t

    sl = soln_t(bp%gr, CMPLX(omega, KIND=WP), c_ext_t(discrim))

    do i = 1, 4
       call sl%set_y(i, y_c(i,:), dy_c_dx(i,:))
    end do

    call sl%set_y(5, SPREAD(CMPLX(0._WP, KIND=WP), 1, bp%n_k), &
                     SPREAD(CMPLX(0._WP, KIND=WP), 1, bp%n_k))

    ! Finish

    return

  end function soln_t_

end module gyre_ad_bvp
