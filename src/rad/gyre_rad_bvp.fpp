! Module   : gyre_rad_bvp
! Purpose  : radial adiabatic bounary value problem solver
!
! Copyright 2013-2019 Rich Townsend
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
  use gyre_interp
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

    if (os_p%gamma_gr /= 1._WP) then
       $WARN(gamma_gr is ignored in 2nd-order radial equations)
    endif
    
    pt_i = gr%pt_i()
    pt_o = gr%pt_o()

    ! Initialize the boundary conditions

    bd = rad_bound_t(cx, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_k-1))

    !$OMP PARALLEL DO
    do k = 1, gr%n_k-1
       df(k) = rad_diff_t(cx, gr%pt(k), gr%pt(k+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t_(bd, df, nm_p) 

    ! Other initializations

    bp%cx => cx
    bp%gr = gr

    bp%tr = rad_trans_t(cx, md_p, os_p)
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

  function wave_t_hom_ (bp, st, j) result (wv)

    class(rad_bvp_t), intent(inout) :: bp
    type(r_state_t), intent(in)     :: st
    integer, intent(in)             :: j
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

    wv = wave_t_y_(bp, st, y, j)

    ! Finish

    return

  end function wave_t_hom_

  !****

  function wave_t_inhom_ (bp, st, z_i, z_o, j) result (wv)

    class(rad_bvp_t), intent(inout) :: bp
    type(r_state_t), intent(in)     :: st
    real(WP), intent(in)            :: z_i(:)
    real(WP), intent(in)            :: z_o(:)
    integer, intent(in)             :: j
    type(wave_t)                    :: wv

    real(WP) :: y(2,bp%n_k)
    integer  :: k

    $CHECK_BOUNDS(SIZE(z_i),bp%n_i)
    $CHECK_BOUNDS(SIZE(z_o),bp%n_o)

    ! Calculate the solution vector

    call bp%build(st)
    call bp%factor()

    y = bp%soln_vec_inhom(z_i, z_o)

    ! Convert to canonical form

    !$OMP PARALLEL DO
    do k = 1, bp%n_k
       call bp%tr%trans_vars(y(:,k), k, st, from=.FALSE.)
    end do

    ! Construct the wave_t

    wv = wave_t_y_(bp, st, y, j)

    ! Finish

    return

  end function wave_t_inhom_

  !****

  function wave_t_y_ (bp, st, y, j) result (wv)

    class(rad_bvp_t), intent(inout) :: bp
    type(r_state_t), intent(in)     :: st
    real(WP), intent(in)            :: y(:,:)
    integer, intent(in)             :: j
    type(wave_t)                    :: wv

    class(model_t), pointer :: ml
    integer                 :: k
    real(WP)                :: U
    real(WP)                :: c_1(bp%n_k)
    real(WP)                :: y_4(bp%n_k)
    real(WP)                :: deul_phi(bp%n_k)
    integer                 :: s
    integer                 :: k_i
    integer                 :: k_o
    type(r_interp_t)        :: in
    real(WP)                :: eul_phi(bp%n_k)
    real(WP)                :: y_3(bp%n_k)
    real(WP)                :: y_g(4,bp%n_k)
    complex(WP)             :: y_c(6,bp%n_k)
    type(c_state_t)         :: st_c
    type(c_ext_t)           :: discrim

    ! Set up gravitational eigenfunctions

    ml => bp%cx%model()

    !$OMP PARALLEL DO PRIVATE (U)
    do k = 1, bp%n_k

       associate ( pt => bp%gr%pt(k) )

         U = ml%coeff(I_U, pt)
         c_1(k) = ml%coeff(I_C_1, pt)

         ! Evaluate y_4
         
         y_4(k) = -U*y(1,k)

         ! Evaluate the Eulerian potential gradient (gravity) perturbation

         if (pt%x /= 0._WP) then
            deul_phi(k) = y_4(k)/(c_1(k)*pt%x)
         else
            deul_phi(k) = 0._WP
         end if

       end associate

    end do

    ! Evaluate the Eulerian potential perturbation,
    ! segment-by-segment, by integrating an interpolant fit to the
    ! gravity perturbation

    seg_loop : do s = bp%gr%s_i(), bp%gr%s_o()

       k_i = bp%gr%k_s_i(s)
       k_o = bp%gr%k_s_o(s)

       in = r_interp_t(bp%gr%pt(k_i:k_o)%x, deul_phi(k_i:k_o), 'MONO')
       
       eul_phi(k_i:k_o) = in%int_f()

    end do seg_loop

    ! Adjust the potential perturbation so that it satisfies the
    ! surface boundary condition
    
    eul_phi = eul_phi - eul_phi(bp%n_k) - y_4(bp%n_k)

    ! Set up y_3 based on it

    y_3 = c_1*eul_phi

    ! Store eigenfunctions into y_g

    y_g(1:2,:) = y

    y_g(3,:) = y_3
    y_g(4,:) = y_4

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

    wv = wave_t(st_c, y_c, discrim, bp%cx, bp%gr, bp%md_p, bp%os_p, j)

    ! Finish

    return

  end function wave_t_y_

end module gyre_rad_bvp
