! Module   : gyre_rad_bvp
! Purpose  : adiabatic radial bounary value problem solver
!
! Copyright 2013-2018 Rich Townsend
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

module gyre_rad_bvp

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_context
  use gyre_ext
  use gyre_grid
  use gyre_model
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_qad_eval
  use gyre_rad_bound
  use gyre_rad_diff
  use gyre_rad_trans
  use gyre_state
  use gyre_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_bvp_t) :: rad_bvp_t
     private
     type(context_t), pointer :: cx => null()
     type(grid_t)             :: gr
     type(rad_trans_t)        :: tr
     type(qad_eval_t)         :: qe
     type(mode_par_t)         :: md_p
     type(osc_par_t)          :: os_p
  end type rad_bvp_t

  ! Interfaces

  interface rad_bvp_t
     module procedure rad_bvp_t_
  end interface rad_bvp_t

  interface wave_t
     module procedure wave_t_hom_
     module procedure wave_t_inhom_
  end interface wave_t

  ! Access specifiers

  private

  public :: rad_bvp_t
  public :: wave_t

  ! Procedures

contains

  function rad_bvp_t_ (cx, gr, md_p, nm_p, os_p) result (bp)

    type(context_t), pointer, intent(in) :: cx
    type(grid_t), intent(in)             :: gr
    type(mode_par_t), intent(in)         :: md_p
    type(num_par_t), intent(in)          :: nm_p
    type(osc_par_t), intent(in)          :: os_p
    type(rad_bvp_t)                      :: bp

    type(point_t)                 :: pt_i
    type(point_t)                 :: pt_o
    type(rad_bound_t)             :: bd
    integer                       :: k
    type(rad_diff_t), allocatable :: df(:)
    type(osc_par_t)               :: qad_os_p

    ! Construct the rad_bvp_t

    if (os_p%cowling_approx) then
       $WARN(cowling_approx is ignored in 2nd-order radial equations)
    endif
    
    pt_i = gr%pt(1)
    pt_o = gr%pt(gr%n_k)

    ! Initialize the boundary conditions

    bd = rad_bound_t(cx, pt_i, pt_o, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_k-1))

    !$OMP PARALLEL DO
    do k = 1, gr%n_k-1
       df(k) = rad_diff_t(cx, pt_i, gr%pt(k), gr%pt(k+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t_(bd, df, nm_p) 

    ! Other initializations

    bp%cx => cx
    bp%gr = gr

    bp%tr = rad_trans_t(cx, pt_i, md_p, os_p)
    call bp%tr%stencil(gr%pt)

    if (os_p%quasiad_eigfuncs) then
       qad_os_p = os_p
       qad_os_p%variables_set = 'GYRE'
       bp%qe = qad_eval_t(cx, gr, md_p, qad_os_p)
    endif

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function rad_bvp_t_

  !****

  function wave_t_hom_ (bp, st) result (wv)

    class(rad_bvp_t), intent(inout) :: bp
    type(r_state_t), intent(in)     :: st
    type(wave_t)                    :: wv

    real(WP) :: y(2,bp%n_k)
    integer  :: k

    ! Calculate the solution vector

    call bp%build(st)
    call bp%factor()

    y = bp%soln_vec_hom()

    ! Convert to canonical form

    !$OMP PARALLEL DO
    do k = 1, bp%n_k
       call bp%tr%trans_vars(y(:,k), k, st, from=.FALSE.)
    end do

    ! Construct the wave_t

    wv = wave_t_y_(bp, st, y)

    ! Finish

    return

  end function wave_t_hom_

  !****

  function wave_t_inhom_ (bp, st, w_i, w_o) result (wv)

    class(rad_bvp_t), intent(inout) :: bp
    type(r_state_t), intent(in)     :: st
    real(WP), intent(in)            :: w_i(:)
    real(WP), intent(in)            :: w_o(:)
    type(wave_t)                    :: wv

    real(WP) :: y(2,bp%n_k)
    integer  :: k

    $CHECK_BOUNDS(SIZE(w_i),bp%n_i)
    $CHECK_BOUNDS(SIZE(w_o),bp%n_o)

    ! Calculate the solution vector

    call bp%build(st)
    call bp%factor()

    y = bp%soln_vec_inhom(w_i, w_o)

    ! Convert to canonical form

    !$OMP PARALLEL DO
    do k = 1, bp%n_k
       call bp%tr%trans_vars(y(:,k), k, st, from=.FALSE.)
    end do

    ! Construct the wave_t

    wv = wave_t_y_(bp, st, y)

    ! Finish

    return

  end function wave_t_inhom_

  !****

  function wave_t_y_ (bp, st, y) result (wv)

    class(rad_bvp_t), intent(inout) :: bp
    type(r_state_t), intent(in)     :: st
    real(WP), intent(in)            :: y(:,:)
    type(wave_t)                    :: wv

    real(WP)        :: y_g(4,bp%n_k)
    integer         :: k
    real(WP)        :: U
    complex(WP)     :: y_c(6,bp%n_k)
    type(c_state_t) :: st_c
    type(c_ext_t)   :: discrim

    ! Set up gravitational eigenfunctions

    y_g(1:2,:) = y

    !$OMP PARALLEL DO PRIVATE (U)
    do k = 1, bp%n_k

       associate (pt => bp%gr%pt(k))
         U = bp%cx%ml%coeff(I_U, pt)
       end associate

       y_g(4,k) = -U*y(1,k)
       
    end do

    y_g(3,:) = 0._WP

    ! Set up complex eigenfunctions

    st_c = c_state_t(CMPLX(st%omega, KIND=WP), st%omega)

    if (bp%os_p%quasiad_eigfuncs) then

       y_c = bp%qe%y_qad(st_c, y_g)

    else

       y_c(1:4,:) = y_g
       y_c(5:6,:) = 0._WP

    endif
 
    ! Construct the wave_t

    discrim = c_ext_t(bp%det())

    wv = wave_t(st_c, y_c, discrim, bp%cx, bp%gr, bp%md_p, bp%os_p)

    ! Finish

    return

  end function wave_t_y_

end module gyre_rad_bvp
