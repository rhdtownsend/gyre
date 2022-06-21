! Module   : gyre_wave
! Purpose  : wave function data
!
! Copyright 2013-2022 Rich Townsend & The GYRE Team
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

module gyre_wave

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_constants
  use gyre_context
  use gyre_ext
  use gyre_freq_context
  use gyre_grid
  use gyre_grid_util
  use gyre_math
  use gyre_model
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_state
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: wave_t
     private
     type(c_state_t)          :: st
     type(context_t)          :: cx
     type(grid_t)             :: gr
     type(mode_par_t), public :: md_p
     type(num_par_t), public  :: nm_p
     type(osc_par_t), public  :: os_p
     type(point_t)            :: pt_i
     type(point_t)            :: pt_o
     complex(WP), allocatable :: y_c(:,:)
     real(WP)                 :: E_scl2
     type(c_ext_t), public    :: discrim
     complex(WP), public      :: scl
     complex(WP), public      :: omega
     complex(WP), public      :: l_i
     integer, public          :: n
     integer, public          :: j_ref
     integer, public          :: l
     integer, public          :: m
     logical, public          :: static
     integer, public          :: id
   contains
     private
     procedure, public :: state
     procedure, public :: context
     procedure, public :: grid
     procedure, public :: freq
     procedure, public :: dfreq_rot
     procedure, public :: x
     procedure, public :: y_i
     procedure, public :: xi_r
     procedure, public :: xi_h
     procedure, public :: eul_phi
     procedure, public :: deul_phi
     procedure, public :: lag_S
     procedure, public :: lag_L
     procedure, public :: eul_P
     procedure, public :: lag_P
     procedure, public :: eul_rho
     procedure, public :: lag_rho
     procedure, public :: eul_T
     procedure, public :: lag_T
     procedure, public :: lambda
     procedure, public :: dE_dx
     procedure, public :: dW_dx
     procedure, public :: dW_eps_dx
     procedure, public :: dzeta_dx
     procedure         :: dzeta_dx_pesnell_
     procedure         :: dzeta_dx_kawaler_
     procedure         :: dzeta_dx_kawaler_grav_
     procedure         :: dzeta_dx_dupret_
     procedure, public :: dzeta_dm
     procedure         :: dzeta_dm_kawaler_
     procedure, public :: dbeta_dx
     procedure, public :: dtau_ss_dx
     procedure, public :: dtau_tr_dx
     procedure, public :: Yt_1
     procedure, public :: Yt_2
     procedure, public :: I_0
     procedure, public :: I_1
     procedure, public :: alpha_0
     procedure, public :: alpha_1
     procedure, public :: prop_type
     procedure, public :: lag_T_eff
     procedure, public :: lag_g_eff
     procedure, public :: f_T
     procedure, public :: f_g
     procedure, public :: psi_T
     procedure, public :: psi_g
     procedure, public :: E
     procedure, public :: E_p
     procedure, public :: E_g
     procedure, public :: E_norm
     procedure, public :: E_ratio
     procedure, public :: H
     procedure, public :: W
     procedure, public :: W_eps
     procedure, public :: tau_ss
     procedure, public :: tau_tr
     procedure, public :: omega_int
     procedure, public :: beta
     procedure, public :: dOmega_rot
     procedure, public :: eta
  end type wave_t

  ! Interfaces

  interface wave_t
     module procedure wave_t_
  end interface wave_t

  interface reallocate
     module procedure reallocate_1_
  end interface reallocate

  ! Access specifiers

  private

  public :: wave_t
  public :: reallocate

  ! Procedures

contains

  function wave_t_ (st, y_c, discrim, cx, gr, md_p, nm_p, os_p, id, static) result (wv)

    type(c_state_t), intent(in)   :: st
    complex(WP), intent(in)       :: y_c(:,:)
    type(c_ext_t), intent(in)     :: discrim
    type(context_t), intent(in)   :: cx
    type(grid_t), intent(in)      :: gr
    type(mode_par_t), intent(in)  :: md_p
    type(num_par_t), intent(in)   :: nm_p
    type(osc_par_t), intent(in)   :: os_p
    integer, intent(in)           :: id
    logical, intent(in), optional :: static
    type(wave_t)                  :: wv

    real(WP) :: x_ref

    $CHECK_BOUNDS(SIZE(y_c, 1),6)
    $CHECK_BOUNDS(SIZE(y_c, 2),gr%n)

    ! Construct the wave_t

    wv%st = st
    wv%cx = cx
    wv%gr = gr

    wv%md_p = md_p
    wv%nm_p = nm_p
    wv%os_p = os_p

    wv%y_c = y_c

    wv%omega = st%omega
    wv%discrim = discrim

    wv%scl = 1._WP

    associate (ml => cx%model(), &
               pt_i => cx%point_i())
      wv%l_i = cx%l_e(cx%Omega_rot(pt_i), st)
    end associate

    wv%n = gr%n

    wv%l = md_p%l
    wv%m = md_p%m
    
    wv%id = id

    if (PRESENT(static)) then
       wv%static = static
    else
       wv%static = .FALSE.
    end if

    ! Locate the reference point

    wv%pt_i = gr%pt_i()
    wv%pt_o = gr%pt_o()

    x_ref = MIN(MAX(os_p%x_ref, wv%pt_i%x), wv%pt_o%x)

    wv%j_ref = MINLOC(abs(gr%pt%x - x_ref), DIM=1)

    ! Cache the ratio E/|scl|**2, to avoid having to perform (expensive)
    ! E integrations in the future

    wv%E_scl2 = wv%E(use_cache=.FALSE.)/abs(wv%scl)**2

    ! Finish

    return

  end function wave_t_

  !****

  $REALLOCATE(type(wave_t),1)

  !****

  function state (this) result (st)

    class(wave_t), intent(in) :: this
    type(c_state_t)           :: st

    ! Return the wave's state

    st = this%st

    ! Finish

    return

  end function state

  !****

  function context (this) result (cx)

    class(wave_t), intent(in) :: this
    type(context_t)           :: cx

    ! Return the wave's context

    cx = this%cx

    ! Finish

    return

  end function context

  !****

  function grid (this) result (gr)

    class(wave_t), intent(in) :: this
    type(grid_t)              :: gr

    ! Return the wave's grid

    gr = this%gr

    ! Finish

    return

  end function grid

  !****

  function freq (this, freq_units, freq_frame)

    class(wave_t), intent(in)          :: this
    character(*), intent(in)           :: freq_units
    character(*), optional, intent(in) :: freq_frame
    complex(WP)                        :: freq

    ! Calculate the frequency

    associate( &
         ml => this%cx%model() )

      if (PRESENT(freq_frame)) then
         freq = (this%st%omega + freq_shift(freq_frame, this%cx, this%md_p))*freq_scale(freq_units, this%cx, this%md_p, this%os_p)
      else
         freq = this%st%omega*freq_scale(freq_units, this%cx, this%md_p, this%os_p)
      endif

    end associate

    ! Finish
    
    return

  end function freq

  !****

  function dfreq_rot (this, freq_units)

    class(wave_t), intent(in) :: this
    character(*), intent(in)  :: freq_units
    real(WP)                  :: dfreq_rot

    associate( &
         ml => this%cx%model() )

      dfreq_rot = this%dOmega_rot()*freq_scale(freq_units, this%cx, this%md_p, this%os_p)
    
    end associate

  end function dfreq_rot

  !****

  function x (this, j)

     class(wave_t), intent(in) :: this
     integer, intent(in)       :: j
     real(WP)                  :: x

     ! Return the abscissa

     x = this%gr%pt(j)%x

     ! Finish

     return

  end function x

  !****

  function y_i (this, i, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: i
    integer, intent(in)       :: j
    complex(WP)               :: y_i

    ! Evaluate y(i)

    y_i = this%scl*this%y_c(i, j)

    ! Finish

    return

  end function y_i

  !****

  function xi_r (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: xi_r

    real(WP)    :: x
    complex(WP) :: y_1

    ! Evaluate the radial displacement perturbation, in units of
    ! R_star

    associate ( &
         l_i => this%l_i )

      x = this%x(j)

      y_1 = this%y_i(1, j)

      if (l_i /= 1._WP) then

         if (x /= 0._WP) then
            xi_r = y_1*pow(x, l_i-1._WP)
         else
            xi_r = 0._WP
         endif
       
      else
       
         xi_r = y_1

      endif

    end associate

    ! Finish

    return

  end function xi_r

  !****

  function xi_h (this, j)
    
    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: xi_h

    complex(WP) :: y_2
    complex(WP) :: y_3
    complex(WP) :: y_4
    real(WP)    :: U
    real(WP)    :: c_1 
    real(WP)    :: Omega_rot
    complex(WP) :: omega_c
    complex(WP) :: lambda

    ! Evaluate the horizontal displacement perturbation, in units of
    ! R_star

    if (this%l /= 0) then

       associate ( &
            ml => this%cx%model(), &
            pt => this%gr%pt(j), &
            l_i => this%l_i)

         y_3 = this%y_i(3, j)

         Omega_rot = this%cx%Omega_rot(pt)

         if (this%static) then

            ! Static cases require a slightly different approach to
            ! evaluate xi_h, since (y_2+y_3)/(c_1*omega_c**2) = 0/0 is
            ! undefined. Instead, use the continuity equation under
            ! the assumption of incompressibility

            y_4 = this%y_i(4, j)

            U = ml%coeff(I_U, pt)

            lambda = this%cx%lambda(Omega_rot, this%st)

            if (l_i /= 1._WP) then

               if (pt%x /= 0._WP) then
                  xi_h = -((4._WP-U)*y_3+y_4)*pow(pt%x, l_i-1._WP)/lambda
               else
                  xi_h = 0._WP
               end if

            else
               
               xi_h = -((4._WP-U)*y_3+y_4)/lambda

            endif

         else

            ! Non-static modes

            y_2 = this%y_i(2, j)

            c_1 = ml%coeff(I_C_1, pt)

            omega_c = this%cx%omega_c(Omega_rot, this%st)

            if (l_i /= 1._WP) then

               if (pt%x /= 0._WP) then
                  xi_h = (y_2+y_3)*pow(pt%x, l_i-1._WP)/(c_1*omega_c**2)
               else
                  xi_h = 0._WP
               end if
               
            else
            
               xi_h = (y_2+y_3)/(c_1*omega_c**2)

            endif

         endif

       end associate

    else

       xi_h = 0._WP

    end if

    ! Finish

    return
    
  end function xi_h

  !****

  function eul_phi (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: eul_phi

    complex(WP) :: y_3
    real(WP)    :: c_1
    
    ! Evaluate the Eulerian gravitational potential perturbation, in
    ! units of G M_star / R_star

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j), &
         l_i => this%l_i )

      y_3 = this%y_i(3, j)

      c_1 = ml%coeff(I_C_1, pt)

      if (l_i /= 0._WP) then

         if (pt%x /= 0._WP) then
            eul_phi = y_3*pow(pt%x, l_i)/c_1
         else
            eul_phi = 0._WP
         endif

      else

         eul_phi = y_3/c_1

      endif

    end associate

    ! Finish

    return

  end function eul_phi

  !****

  function deul_phi (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: deul_phi

    complex(WP) :: y_4
    real(WP)    :: c_1
    
    ! Evaluate the Eulerian potential gradient (gravity) perturbation,
    ! in units of G M_star / R_star**2

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j), &
         l_i => this%l_i )

      y_4 = this%y_i(4, j)

      c_1 = ml%coeff(I_C_1, pt)

      if (l_i /= 1._WP) then

         if (pt%x /= 0._WP) then
            deul_phi = y_4*pow(pt%x, l_i-1._WP)/c_1
         else
            deul_phi = 0._WP
         end if
       
      else

         deul_phi = y_4/c_1

      end if

    end associate

    ! Finish

    return

  end function deul_phi

  !****

  function lag_S (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: lag_S

    complex(WP) :: y_5

    ! Evaluate the Lagrangian specific entropy perturbation, in units
    ! of c_p

    associate ( &
         pt => this%gr%pt(j), &
         l_i => this%l_i )

      y_5 = this%y_i(5, j)

      if (pt%x /= 0._WP) then
         lag_S = y_5*pow(pt%x, l_i-2._WP)
      else
         lag_S = 0._WP
      endif

    end associate

    ! Finish

    return

  end function lag_S

  !****

  function lag_L (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: lag_L

    complex(WP) :: y_6

    ! Evaluate the Lagrangian radiative luminosity perturbation, in
    ! units of L_star

    associate ( &
         pt => this%gr%pt(j), &
         l_i => this%l_i )

      y_6 = this%y_i(6, j)

      if (pt%x /= 0._WP) then
         lag_L = y_6*pow(pt%x, l_i+1._WP)
      else
         lag_L = 0._WP
      endif

    end associate

    ! Finish

    return

  end function lag_L

  !****

  function eul_P (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: eul_P

    complex(WP) :: xi_r
    complex(WP) :: lag_P
    real(WP)    :: V_2

    ! Evaluate the Eulerian pressure perturbation, in units of P

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      xi_r = this%xi_r(j)
      lag_P = this%lag_P(j)

      V_2 = ml%coeff(I_V_2, pt)

      eul_P = lag_P + V_2*pt%x*xi_r

    end associate

    ! Finish

    return

  end function eul_P
  
  !****

  function lag_P (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: lag_P

    complex(WP) :: y_1
    complex(WP) :: y_2
    real(WP)    :: V_2

    ! Evaluate the Lagrangian pressure perturbation, in units of P

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j), &
         l_i => this%l_i )

      y_1 = this%y_i(1, j)
      y_2 = this%y_i(2, j)

      V_2 = ml%coeff(I_V_2, pt)

      if (l_i /= 0._WP) then

         if (pt%x /= 0._WP) then
            lag_P = V_2*(y_2 - y_1)*pow(pt%x, l_i)
         else
            lag_P = 0._WP
         end if

      else

         lag_P = V_2*(y_2 - y_1)

      endif

    end associate

    ! Finish

    return

  end function lag_P

  !****

  function eul_rho (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: eul_rho

    complex(WP) :: xi_r
    complex(WP) :: lag_rho
    real(WP)    :: U
    real(WP)    :: dU
    real(WP)    :: D

    ! Evaluate the Eulerian density perturbation, in units of rho

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      xi_r = this%xi_r(j)
      lag_rho = this%lag_rho(j)

      U = ml%coeff(I_U, pt)
      dU =ml%dcoeff(I_U, pt)

      D = dU + U - 3._WP

      if (pt%x /= 0._WP) then
         eul_rho = lag_rho - D*xi_r/pt%x
      else
         eul_rho = lag_rho
      endif

    end associate

    ! Finish

    return

  end function eul_rho

  !****

  function lag_rho (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: lag_rho

    complex(WP) :: lag_P
    complex(WP) :: lag_S
    real(WP)    :: Gamma_1
    real(WP)    :: delta

    ! Evaluate the Lagrangian density perturbation, in units of
    ! rho. This expression implements eqn. (13.83) of [Unno:1989]

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      lag_P = this%lag_P(j)
      lag_S = this%lag_S(j)

      Gamma_1 = ml%coeff(I_GAMMA_1, pt)
      delta = ml%coeff(I_DELTA, pt)

      lag_rho = lag_P/Gamma_1 - delta*lag_S

    end associate

    ! Finish

    return
    
  end function lag_rho

  !****

  function eul_T (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: eul_T

    complex(WP) :: xi_r
    complex(WP) :: lag_T
    real(WP)    :: V_2
    real(WP)    :: nabla

    ! Evaluate the Lagrangian temperature perturbation, in units of T

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      xi_r = this%xi_r(j)
      lag_T = this%lag_T(j)

      V_2 = ml%coeff(I_V_2, pt)
      nabla = ml%coeff(I_NABLA, pt)
      
      eul_T = lag_T + nabla*V_2*pt%x*xi_r

    end associate

    ! Finish

    return

  end function eul_T

  !****

  function lag_T (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: lag_T

    complex(WP) :: lag_P
    complex(WP) :: lag_S
    real(WP)    :: nabla_ad

    ! Evaluate the Lagrangian temperature perturbation, in units of
    ! T. This expression implements eqn. (13.84) of [Unno:1989]

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      lag_P = this%lag_P(j)
      lag_S = this%lag_S(j)

      nabla_ad = ml%coeff(I_NABLA_AD, pt)
      
      lag_T = nabla_ad*lag_P + lag_S

    end associate

    ! Finish

    return

  end function lag_T

  !****

  function lambda (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: lambda

    real(WP) :: Omega_rot

    ! Evaluate the angular eigenvalue

    associate ( &
         pt => this%gr%pt(j) )
      
      Omega_rot = this%cx%Omega_rot(pt)
    
      lambda = this%cx%lambda(Omega_rot, this%st)

    end associate

    ! Finish

    return

  end function lambda
    
  !****

  function dE_dx (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    real(WP)                  :: dE_dx

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    complex(WP) :: lambda
    real(WP)    :: U
    real(WP)    :: c_1

    ! Evaluate the differential inertia, in units of M_star
    ! R_star**2. This expression is based on eqn. (3.139) of [Aerts:2010]

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      xi_r = this%xi_r(j)
      xi_h = this%xi_h(j)

      lambda = this%lambda(j)

      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      dE_dx = (abs(xi_r)**2 + abs(lambda)*abs(xi_h)**2)*U*pt%x**2/c_1

    end associate

    ! Finish

    return

  end function dE_dx

  !****

  function dW_dx (this, j)

    use gyre_evol_model

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    real(WP)                  :: dW_dx

    real(WP)    :: t_dyn
    real(WP)    :: t_kh
    complex(WP) :: lag_T
    complex(WP) :: lag_S
    real(WP)    :: c_thk

    ! Evaluate the differential work, in units of G M_star**2/R_star.
    ! This expression is based on eqn. (25.9) of [Unno:1989]

    select type (ml => this%cx%model())
    class is (evol_model_t)
       t_dyn = sqrt(ml%R_star**3/(G_GRAVITY*ml%M_star))
       t_kh = (G_GRAVITY*ml%M_star**2/ml%R_star)/ml%L_star
    class default
       t_dyn = 1._WP
       t_kh = 1._WP
    end select

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      lag_T = this%lag_T(j)
      lag_S = this%lag_S(j)
    
      c_thk = ml%coeff(I_C_THK, pt)

      dW_dx = PI*AIMAG(CONJG(lag_T)*lag_S)*c_thk*pt%x**2*t_dyn/t_kh

    end associate

    ! Finish

    return

  end function dW_dx

  !****

  function dW_eps_dx (this, j)

    use gyre_evol_model

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    real(WP)                  :: dW_eps_dx

    real(WP)    :: t_dyn
    real(WP)    :: t_kh
    complex(WP) :: lag_rho
    complex(WP) :: lag_T
    real(WP)    :: c_eps
    complex(WP) :: eps_rho
    complex(WP) :: eps_T
    real(WP)    :: omega_r

    ! Evaluate the differential work associated with nuclear
    ! processes, in units of G M_star**2/R_star.  This expression is
    ! based on eqn. (25.9) of [Unno:1989]
    
    select type (ml => this%cx%model())
    class is (evol_model_t)
       t_dyn = sqrt(ml%R_star**3/(G_GRAVITY*ml%M_star))
       t_kh = (G_GRAVITY*ml%M_star**2/ml%R_star)/ml%L_star
    class default
       t_dyn = 1._WP
       t_kh = 1._WP
    end select

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      lag_rho = this%lag_rho(j)
      lag_T = this%lag_T(j)

      c_eps = ml%coeff(I_C_EPS, pt)
      
      eps_rho = this%cx%eps_rho(this%st, pt)
      eps_T = this%cx%eps_T(this%st, pt)
      
      omega_r = REAL(this%st%omega)

      dW_eps_dx = PI/omega_r*REAL(CONJG(lag_T)*c_eps*(eps_rho*lag_rho + eps_T*lag_T))*pt%x**2*t_dyn/t_kh

    end associate

    ! Finish

    return

  end function dW_eps_dx

  !****

  function dzeta_dx (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: dzeta_dx

    ! Calculate the dimensionless frequency weight function.
    
    select case (this%os_p%zeta_scheme)
    case ('PESNELL')
       dzeta_dx = this%dzeta_dx_pesnell_(j)
    case ('KAWALER_GRAV')
       dzeta_dx = this%dzeta_dx_kawaler_grav_(j)
    case ('KAWALER')
       dzeta_dx = this%dzeta_dx_kawaler_(j)
    case ('DUPRET')
       dzeta_dx = this%dzeta_dx_dupret_(j)
    case default
       $ABORT(Invalid zeta_scheme)
    end select

    ! Finish

    return

  end function dzeta_dx
    
  !****

  function dzeta_dx_pesnell_ (this, j) result (dzeta_dx)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: dzeta_dx

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    complex(WP) :: eul_phi
    complex(WP) :: eul_rho
    complex(WP) :: lag_rho
    complex(WP) :: lag_P
    complex(WP) :: lambda
    real(WP)    :: V_2
    real(WP)    :: U
    real(WP)    :: c_1
    real(WP)    :: x4_V

    ! Calculate the dimensionless frequency weight function.  This is
    ! based on the derivative of equation (A5) of [Pesnell:1987] with
    ! respect to x. NOTE: This doesn't seem to agree in any way with
    ! the other weight functions (Kawaler, Dupret)

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      xi_r = this%xi_r(j)
      xi_h = this%xi_h(j)
      eul_phi = this%eul_phi(j)
      eul_rho = this%eul_rho(j)
      lag_rho = this%lag_rho(j)
      lag_P = this%lag_P(j)

      lambda = this%lambda(j)

      V_2 = ml%coeff(I_V_2, pt)
      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      x4_V = pt%x**2/V_2

      dzeta_dx = CONJG(lag_rho)*lag_P*(U*x4_V/(c_1**2)) + &
                 2._WP*REAL(lambda*CONJG(xi_r)*xi_h/c_1*(pt%x**2*U/c_1)) + &
                 CONJG(eul_rho)*eul_phi*(pt%x**2*U/c_1) + &
                 abs(xi_r)**2*(U-4._WP)/c_1*(pt%x**2*U/c_1)

    end associate

    ! Finish

    return

  end function dzeta_dx_pesnell_

  !****

  function dzeta_dx_kawaler_ (this, j) result (dzeta_dx)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: dzeta_dx

    complex(WP) :: xi_r
    complex(WP) :: eul_P
    complex(WP) :: eul_phi
    complex(WP) :: deul_phi
    complex(WP) :: lambda
    real(WP)    :: V_2
    real(WP)    :: U
    real(WP)    :: As
    real(WP)    :: c_1
    real(WP)    :: Gamma_1
    real(WP)    :: x4_V

    ! Calculate the dimensionless frequency weight function.  This is
    ! based on the derivative of equation (7) of [Kawaler:1985] with
    ! respect to x; note that that equation was derived for adiabatic
    ! pulsation, and the simple extension to non-adiabatic pulsation
    ! implemented here may not be quite correct

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      xi_r = this%xi_r(j)
      eul_P = this%eul_P(j)
      eul_phi = this%eul_phi(j)
      deul_phi = this%deul_phi(j)

      V_2 = ml%coeff(I_V_2, pt)
      As = ml%coeff(I_AS, pt)
      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)
      Gamma_1 = ml%coeff(I_GAMMA_1, pt)
      
      lambda = this%lambda(j)

      x4_V = pt%x**2/V_2

      dzeta_dx = CONJG(eul_P)*eul_P*(U*x4_V/(Gamma_1*c_1**2)) + &
                 CONJG(xi_r)*xi_r*(pt%x**2*U*As/c_1**2) - &
                 CONJG(pt%x*deul_phi + lambda*eul_phi)*(pt%x*deul_phi + lambda*eul_phi)
           
    end associate

    ! Finish

    return

  end function dzeta_dx_kawaler_

  !****

  function dzeta_dx_kawaler_grav_ (this, j) result (dzeta_dx)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: dzeta_dx

    complex(WP) :: xi_r
    complex(WP) :: eul_P
    complex(WP) :: eul_phi
    complex(WP) :: deul_phi
    complex(WP) :: lambda
    real(WP)    :: V_2
    real(WP)    :: U
    real(WP)    :: As
    real(WP)    :: c_1
    real(WP)    :: Gamma_1
    real(WP)    :: x4_V

    ! Calculate the gravitational part of the dimensionless frequency
    ! weight function.  This is based on the derivative of the N term
    ! in equation (7) of [Kawaler:1985] with respect to x; note that
    ! that equation was derived for adiabatic pulsation, and the
    ! simple extension to non-adiabatic pulsation implemented here may
    ! not be quite correct

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      xi_r = this%xi_r(j)
      eul_P = this%eul_P(j)
      eul_phi = this%eul_phi(j)
      deul_phi = this%deul_phi(j)

      V_2 = ml%coeff(I_V_2, pt)
      As = ml%coeff(I_AS, pt)
      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)
      Gamma_1 = ml%coeff(I_GAMMA_1, pt)
      
      lambda = this%lambda(j)

      x4_V = pt%x**2/V_2

      dzeta_dx = CONJG(xi_r)*xi_r*(pt%x**2*U*As/c_1**2)
           
    end associate

    ! Finish

    return

  end function dzeta_dx_kawaler_grav_

  !****

  function dzeta_dx_dupret_ (this, j) result (dzeta_dx)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: dzeta_dx

    complex(WP) :: xi_r
    complex(WP) :: eul_phi
    complex(WP) :: eul_rho
    complex(WP) :: lag_rho
    complex(WP) :: lag_P
    real(WP)    :: V_2
    real(WP)    :: V
    real(WP)    :: U
    real(WP)    :: As
    real(WP)    :: c_1
    real(WP)    :: Gamma_1
    real(WP)    :: x4_V

    ! Calculate the dimensionless frequency weight function.  This is
    ! based on the derivative of equation (1.71) of [Dupret:2002a]
    ! with respect to x

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      xi_r = this%xi_r(j)
      eul_phi = this%eul_phi(j)
      eul_rho = this%eul_rho(j)
      lag_rho = this%lag_rho(j)
      lag_P = this%lag_P(j)

      V_2 = ml%coeff(I_V_2, pt)
      V = V_2*pt%x**2
      As = ml%coeff(I_AS, pt)
      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)
      Gamma_1 = ml%coeff(I_GAMMA_1, pt)

      x4_V = pt%x**2/V_2

      dzeta_dx = CONJG(lag_rho)*lag_P*(U*x4_V/(c_1**2)) + &
                 2._WP*REAL(lag_rho*CONJG(xi_r)*(pt%x/c_1)*(pt%x**2*U/c_1)) + &
                 CONJG(eul_rho)*eul_phi*(pt%x**2*U/c_1) - &
                 abs(xi_r)**2*(-V/Gamma_1-As)/c_1*(pt%x**2*U/c_1)

    end associate

    ! Finish

    return

  end function dzeta_dx_dupret_
  
  !****

  function dzeta_dm (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: dzeta_dm

    ! Calculate the dimensionless frequency weight function.
    
    select case (this%os_p%zeta_scheme)
    case ('KAWALER')
       dzeta_dm = this%dzeta_dm_kawaler_(j)
    case default
       $ABORT(Invalid zeta_scheme)
    end select

    ! Finish

    return

  end function dzeta_dm
    
  !****

  function dzeta_dm_kawaler_ (this, j) result (dzeta_dm)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: dzeta_dm

    complex(WP) :: xi_r
    complex(WP) :: eul_P
    complex(WP) :: eul_phi
    complex(WP) :: deul_phi
    complex(WP) :: lambda
    real(WP)    :: V_2
    real(WP)    :: U
    real(WP)    :: As
    real(WP)    :: c_1
    real(WP)    :: Gamma_1

    ! Calculate the dimensionless frequency weight function.  This is
    ! based on the derivative of equation (7) of [Kawaler:1985] with
    ! respect to m; note that that equation was derived for adiabatic
    ! pulsation, and the simple extension to non-adiabatic pulsation
    ! implemented here may not be quite correct

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      xi_r = this%xi_r(j)
      eul_P = this%eul_P(j)
      eul_phi = this%eul_phi(j)
      deul_phi = this%deul_phi(j)

      V_2 = ml%coeff(I_V_2, pt)
      As = ml%coeff(I_AS, pt)
      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)
      Gamma_1 = ml%coeff(I_GAMMA_1, pt)
      
      lambda = this%lambda(j)

      if (pt%x /= 0) then
         dzeta_dm = CONJG(eul_P)*eul_P/(V_2*Gamma_1*c_1) + &
                    CONJG(xi_r)*xi_r*(As/c_1) - &
                    CONJG(pt%x*deul_phi + lambda*eul_phi)*(pt%x*deul_phi + lambda*eul_phi)*(c_1/(U*pt%x**2))
      else
         dzeta_dm = 0._WP
      endif
           
    end associate

    ! Finish

    return

  end function dzeta_dm_kawaler_

  !****

  function dbeta_dx (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    real(WP)                  :: dbeta_dx

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    real(WP)    :: U
    real(WP)    :: c_1
    complex(WP) :: lambda
    real(WP)    :: E
    
    ! Calculate the (unnormalized) rotation splitting kernel. This is
    ! based on the derivative of equation 3.357 of [Aerts:2010] with
    ! respect to x

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      xi_r = this%xi_r(j)
      xi_h = this%xi_h(j)

      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      ! Question: should the following be lambda or l(l+1)?

      lambda = this%lambda(j)

      E = this%E()

      dbeta_dx = REAL((abs(xi_r)**2 + (lambda-1._WP)*abs(xi_h)**2 - &
                       2._WP*xi_r*CONJG(xi_h))*U*pt%x**2/c_1)/E

    end associate

    ! Finish

    return

  end function dbeta_dx

  !****

  function dtau_ss_dx (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    real(WP)                  :: dtau_ss_dx

    complex(WP) :: lag_P
    complex(WP) :: lag_rho
    real(WP)    :: V_2
    real(WP)    :: c_1
    real(WP)    :: U
    
    ! Evaluate the steady-state differential torque, in units of G
    ! M_star**2/R_star. This expression is based on eqn. (13) of
    ! [Townsend:2018]

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j), &
         m => this%m )

      lag_P = this%lag_P(j)

      lag_rho = this%lag_rho(j)

      V_2 = ml%coeff(I_V_2, pt)
      c_1 = ml%coeff(I_C_1, pt)
      U = ml%coeff(I_U, pt)

      dtau_ss_dx = m*pt%x**2*AIMAG(lag_rho*CONJG(lag_P))*(U/(2._WP*c_1**2*V_2))
      
    end associate

    ! Finish

    return

  end function dtau_ss_dx

  !****

  function dtau_tr_dx (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    real(WP)                  :: dtau_tr_dx

    complex(WP) :: xi_r
    complex(WP) :: eul_P
    complex(WP) :: eul_rho
    complex(WP) :: lag_rho
    complex(WP) :: eul_phi
    real(WP)    :: V_2
    real(WP)    :: U
    real(WP)    :: c_1
    real(WP)    :: Omega_rot
    complex(WP) :: omega_c

    ! Evaluate the transient differential torque, in units of G
    ! M_star**2/R_star. This expression is based on eqn. (14) of
    ! [Towsend:2018]

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j), &
         m => this%m )

      xi_r = this%xi_r(j)

      eul_P = this%eul_P(j)

      lag_rho = this%lag_rho(j)
      eul_rho = this%eul_rho(j)

      eul_phi = this%eul_phi(j)

      V_2 = ml%coeff(I_V_2, pt)
      c_1 = ml%coeff(I_C_1, pt)
      U = ml%coeff(I_U, pt)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, this%st)

      dtau_tr_dx = m*pt%x**2*AIMAG((omega_c/CONJG(omega_c) - 1._WP)*( &
           lag_rho*CONJG(eul_P)/(c_1*V_2) + &
           eul_rho*CONJG(eul_phi) + &
           xi_r*CONJG(eul_rho)*pt%x/c_1))*(U/(2._WP*c_1))
           
    end associate

    ! Finish

    return
    
  end function dtau_tr_dx

  !****

  function Yt_1 (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: Yt_1

    complex(WP) :: y_1
    complex(WP) :: y_3
    complex(WP) :: y_4
    real(WP)    :: Jt

    ! Evaluate the Takata Y_1 function. This expression is equivalent to
    ! eqn. (69) of [Takata:2006b], divided by x**(2-l)

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      y_1 = this%y_i(1, j)
      y_3 = this%y_i(3, j)
      y_4 = this%y_i(4, j)

      Jt = 1._WP - ml%coeff(I_U, pt)/3._WP

      Yt_1 = Jt*y_1 + (y_3 - y_4)/3._WP

    end associate

    ! Finish

    return

  end function Yt_1

  !****

  function Yt_2 (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: Yt_2

    complex(WP) :: y_1
    complex(WP) :: y_2

    ! Evaluate the Takata Y_2 function. This expression is equivalent to 
    ! eqn. (70) of [Takata:2006b], divided by V

    y_1 = this%y_i(1, j)
    y_2 = this%y_i(2, j)

    Yt_2 = y_2 - y_1

    ! Finish

    return

  end function Yt_2

  !****

  function I_0 (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: I_0

    complex(WP) :: y_1
    complex(WP) :: y_4
    real(WP)    :: U
    real(WP)    :: c_1
    
    ! Evaluate the I_0 integral, which should be zero for radial
    ! cases. This expression is based on eqn. (42) of [Takata:2006a]

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j), &
         l_i => this%l_i )

      y_1 = this%y_i(1, j)
      y_4 = this%y_i(4, j)

      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      if (pt%x /= 0._WP) then
         I_0 = pow(pt%x, l_i+1._WP)*(U*y_1 + y_4)/c_1
      else
         I_0 = 0._WP
      endif

    end associate

    ! Finish

    return

  end function I_0

  !****

  function I_1 (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: I_1

    complex(WP) :: y_1
    complex(WP) :: y_2
    complex(WP) :: y_3
    complex(WP) :: y_4
    real(WP)    :: U
    real(WP)    :: c_1
    real(WP)    :: Omega_rot
    complex(WP) :: omega_c

    ! Evaluate the I_0 integral, which should be zero for dipole
    ! cases. This expression is based on eqn. (43) of [Takata:2006a]

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j), &
         l_i => this%l_i )

      y_1 = this%y_i(1, j)
      y_2 = this%y_i(2, j)
      y_3 = this%y_i(3, j)
      y_4 = this%y_i(4, j)

      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, this%st)

      if (pt%x /= 0._WP) then
         I_1 = pow(pt%x, l_i+2._WP)*(c_1*omega_c**2*U*y_1 - U*y_2 - &
               (c_1*omega_c**2 + 2._WP)*y_3 + (c_1*omega_c**2 - 1._WP)*y_4)/c_1**2
      else
         I_1 = 0._WP
      endif

    end associate
      
    ! Finish

    return

  end function I_1

  !****

  function alpha_0 (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    real(WP)                  :: alpha_0

    real(WP)    :: U
    real(WP)    :: c_1
    real(WP)    :: Omega_rot
    complex(WP) :: omega_c
    complex(WP) :: lambda

    ! Evaluate the alpha_0 excitation parameter defined in equation
    ! (26.10) of Unno et al. (2017)

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j))

      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, this%st)
      lambda = this%cx%lambda(Omega_rot, this%st)

      alpha_0 = 4._WP - U - abs(lambda)/(c_1*abs(omega_c)**2) + c_1*abs(omega_c)**2

    end associate

    ! Finish

    return

  end function alpha_0

  !****

  function alpha_1 (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    real(WP)                  :: alpha_1

    real(WP)    :: alpha_0
    real(WP)    :: V
    real(WP)    :: U
    real(WP)    :: c_1
    real(WP)    :: Gamma_1
    real(WP)    :: nabla_ad
    real(WP)    :: kap_T
    real(WP)    :: kap_rho
    real(WP)    :: Omega_rot
    complex(WP) :: omega_c
    complex(WP) :: lambda

    ! Evaluate the alpha_1 excitation parameter defined in equation
    ! (26.12) of Unno et al. (2017)

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j))

      alpha_0 = this%alpha_0(j)

      V = ml%coeff(I_V_2, pt)*pt%x**2
      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      Gamma_1 = ml%coeff(I_GAMMA_1, pt)
      nabla_ad = ml%coeff(I_NABLA_AD, pt)

      kap_T = ml%coeff(I_KAP_T, pt)
      kap_rho = ml%coeff(I_KAP_RHO, pt)

      Omega_rot = this%cx%Omega_rot(pt)

      omega_c = this%cx%omega_c(Omega_rot, this%st)
      lambda = this%cx%lambda(Omega_rot, this%st)

      if (V /= 0._WP) then
         alpha_1 = 4._WP - 1._WP/nabla_ad - kap_T - kap_rho/(nabla_ad*Gamma_1) + &
                   (1._WP - abs(lambda)/(c_1*abs(omega_c)**2*V))*(c_1*abs(omega_c)**2 - U)/(nabla_ad*alpha_0)
      else
         alpha_1 = SIGN(HUGE(0._WP), (c_1*abs(omega_c)**2 - U)/(nabla_ad*alpha_0))
      endif
         
    end associate

    ! Finish

    return

  end function alpha_1

  !****

  function prop_type (this, j)

    class(wave_t), intent(in) :: this
    integer, intent(in)       :: j
    integer                   :: prop_type

    real(WP) :: V
    real(WP) :: As
    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: Gamma_1
    real(WP) :: Omega_rot
    real(WP) :: lambda
    real(WP) :: omega_c
    real(WP) :: g_4
    real(WP) :: g_2
    real(WP) :: g_0
    real(WP) :: gamma

    ! Set up the propagation type (0 -> evanescent, 1 -> p, -1 -> g)

    associate ( &
         ml => this%cx%model(), &
         pt => this%gr%pt(j) )

      if (ml%is_vacuum(pt)) then

         prop_type = 0

      else

         ! Calculate the discriminant gamma

         V = ml%coeff(I_V_2, pt)*pt%x**2
         As = ml%coeff(I_AS, pt)
         U = ml%coeff(I_U, pt)
         c_1 = ml%coeff(I_C_1, pt)
         Gamma_1 = ml%coeff(I_GAMMA_1, pt)

         Omega_rot = this%cx%Omega_rot(pt)

         lambda = REAL(this%lambda(j))

         omega_c = REAL(this%cx%omega_c(Omega_rot, this%st))

         g_4 = -4._WP*V/Gamma_1*c_1
         g_2 = (As - V/Gamma_1 - U + 4._WP)**2 + 4._WP*V/Gamma_1*As + 4._WP*lambda
         g_0 = -4._WP*lambda*As/c_1
       
         gamma = (g_4*omega_c**4 + g_2*omega_c**2 + g_0)/omega_c**2

         ! Use the sign of gamma to set up prop_type
         
         if (gamma > 0._WP) then

            prop_type = 0

         else

            if (REAL(c_1*omega_c**2) < As) then
               prop_type = -1
            else
               prop_type = 1
            endif

         endif

      endif

    end associate
    
    ! Finish

    return

  end function prop_type

  !****

  function lag_T_eff (this)

    class(wave_t), intent(in) :: this
    complex(WP)               :: lag_T_eff

    complex(WP) :: xi_r
    complex(WP) :: lag_L

    ! Calculate the effective temperature perturbation at x_ref
    ! (assumed to correspond to the photosphere), in units of
    ! T_eff. This expression is based on the standard definition of
    ! effective temperature

    associate (j => this%j_ref)

      xi_r = this%xi_r(j)
      lag_L = this%lag_L(j)

      lag_T_eff = 0.25_WP*(lag_L - 2._WP*xi_r)

    end associate

    ! Finish

    return

  end function lag_T_eff

  !****

  function lag_g_eff (this)

    class(wave_t), intent(in) :: this
    complex(WP)               :: lag_g_eff

    complex(WP) :: xi_r
    complex(WP) :: deul_phi
    real(WP)    :: c_1
    real(WP)    :: U

    ! Calculate the effective gravity perturbation at x_ref (assumed
    ! to correspond to the photosphere), in units of the gravity. This
    ! expression is based on eqn. (24) of [Dupret:2002b]

    associate ( &
         j => this%j_ref, &
         ml => this%cx%model(), &
         pt => this%gr%pt(this%j_ref), &
         omega => this%st%omega )

      xi_r = this%xi_r(j)
      deul_phi = this%deul_phi(j)

      c_1 = ml%coeff(I_C_1, pt)
      U = ml%coeff(I_U, pt)

      lag_g_eff = (c_1/pt%x)*deul_phi + (U - (2._WP + c_1*omega**2))*xi_r/pt%x

    end associate

    ! Finish

    return

  end function lag_g_eff

  !****

  function f_T (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: f_T

    complex(WP) :: C_T

    ! Evaluate the non-adiabatic f_T parameter. This is expression is
    ! based on eqn. (5) of [Dupret:2003]

    associate (j => this%j_ref)

      C_T = this%lag_T_eff()/this%xi_r(j)

      f_T = abs(C_T)

    end associate

    ! Finish

    return

  end function f_T

  !****

  function f_g (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: f_g

    complex(WP) :: C_g

    ! Evaluate the non-adiabatic f_g parameter. This is expression is
    ! based on eqn. (6) of [Dupret:2003]

    associate (j => this%j_ref)

      C_g = this%lag_g_eff()/this%xi_r(j)

      f_g = -abs(C_g)

    end associate

    ! Finish

    return

  end function f_g

  !****

  function psi_T (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: psi_T

    complex(WP) :: C_T

    ! Calculate the non-adiabatic psi_T parameter, in radians. This is
    ! expression is based on eqn. (5) of [Dup2003]

    associate (j => this%j_ref)

      C_T = this%lag_T_eff()/this%xi_r(j)

      psi_T = atan2(AIMAG(C_T), REAL(C_T))

    end associate

    ! Finish

    return

  end function psi_T

  !****

  function psi_g (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: psi_g

    ! Calculate the non-adiabatic psi_g parameter, in radians

    psi_g = PI

    ! Finish

    return

  end function psi_g

  !****

  function E (this, use_cache)

    class(wave_t), intent(in)     :: this
    logical, intent(in), optional :: use_cache
    real(WP)                      :: E

    logical  :: use_cache_
    integer  :: j
    real(WP) :: x(this%n)
    real(WP) :: dE_dx(this%n)

    ! Calculate the inertia, in units of M_star R_star**2

    if (PRESENT(use_cache)) then
       use_cache_ = use_cache
    else
       use_cache_ = .TRUE.
    endif

    if (use_cache_) then

       E = this%E_scl2*abs(this%scl)**2

    else

       !$OMP PARALLEL DO
       do j = 1, this%n
          x(j) = this%gr%pt(j)%x
          dE_dx(j) = this%dE_dx(j)
       end do

       E = integrate(x, dE_dx)

    end if

    ! Finish

    return

  end function E

  !****

  function E_p (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: E_p

    integer  :: j
    real(WP) :: dE_dx(this%n)

    ! Calculate the inertia in acoustic-wave propagation regions, in
    ! units of M_star R_star**2

    !$OMP PARALLEL DO
    do j = 1, this%n
       if (this%prop_type(j) == 1) then
          dE_dx(j) = this%dE_dx(j)
       else
          dE_dx(j) = 0._WP
       endif
    end do

    E_p = integrate(this%gr%pt%x, dE_dx)

    ! Finish

    return

  end function E_p

  !****

  function E_g (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: E_g

    integer  :: j
    real(WP) :: dE_dx(this%n)

    ! Calculate the inertia in gravity-wave propagation regions, in
    ! units of M_star R_star**2

    !$OMP PARALLEL DO
    do j = 1, this%n
       if (this%prop_type(j) == -1) then
          dE_dx(j) = this%dE_dx(j)
       else
          dE_dx(j) = 0._WP
       endif
    end do

    E_g = integrate(this%gr%pt%x, dE_dx)

    ! Finish

    return

  end function E_g

  !****

  function E_norm (this)
 
    class(wave_t), intent(in) :: this
    real(WP)                  :: E_norm

    real(WP)    :: E
    complex(WP) :: xi_r
    complex(WP) :: xi_h
    complex(WP) :: lambda
    real(WP)    :: A2

    ! Calculate the normalized inertia. This expression is based on
    ! eqn. (3.140) of [Aer2010]

    associate (j => this%j_ref)

      E = this%E()

      xi_r = this%xi_r(j)
      xi_h = this%xi_h(j)

      lambda = this%lambda(j)

      select case (this%os_p%inertia_norm)
      case ('RADIAL')
         A2 = abs(xi_r)**2
      case ('HORIZ')
         A2 = abs(lambda)*abs(xi_h)**2
      case ('BOTH')
         A2 = abs(xi_r)**2 + abs(lambda)*abs(xi_h)**2
      case default
         $ABORT(Invalid inertia_norm)
      end select

      if (A2 == 0._WP) then
         $WARN(Amplitude at x_ref is zero; not normalizing inertia)
         E_norm = E
      else
         E_norm = E/A2
      endif

    end associate

    ! Finish

    return

  end function E_norm

  !****

  function E_ratio (this)
 
    class(wave_t), intent(in) :: this
    real(WP)                  :: E_ratio

    integer  :: j
    real(WP) :: dE_dx(this%n)
    real(WP) :: E_above
    real(WP) :: E_below

    ! Calculate the ratio of the inertias above and below the
    ! reference point

    !$OMP PARALLEL DO
    do j = 1, this%n
       dE_dx(j) = this%dE_dx(j)
    end do

    if (this%j_ref < this%n) then
       associate (j_a => this%j_ref, j_b => this%n)
         E_above = integrate(this%gr%pt(j_a:j_b)%x, dE_dx(j_a:j_b))
       end associate
    else
       E_above = 0._WP
    endif

    if (this%j_ref > 1) then
       associate (j_a => 1, j_b => this%j_ref)
         E_below = integrate(this%gr%pt(j_a:j_b)%x, dE_dx(j_a:j_b))
       end associate
    else
       E_below = 0._WP
    endif

    E_ratio = E_above/(E_below + E_above)

    ! Finish

    return

  end function E_ratio

  !****

  function H (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: H

    ! Calculate the total energy (i.e., twice the volume-integrated
    ! kinetic energy), in units of G M_star**2 / R_star

    H = 0.5_WP*REAL(this%st%omega)**2*this%E()

    ! Finish

    return

  end function H

  !****

  function W (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: W

    integer  :: j
    real(WP) :: dW_dx(this%n)
    
    ! Calculate the total work, in units of G M_star**2/R_star

    !$OMP PARALLEL DO
    do j = 1, this%n
       dW_dx(j) = this%dW_dx(j)
    end do

    W = integrate(this%gr%pt%x, dW_dx)

    ! Finish

    return

  end function W

  !****

  function W_eps (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: W_eps
    
    integer  :: j
    real(WP) :: dW_eps_dx(this%n)
    
    ! Calculate the total work associated with nuclear processes, in
    ! units of G M_star**2/R_star

    !$OMP PARALLEL DO
    do j = 1, this%n
       dW_eps_dx(j) = this%dW_eps_dx(j)
    end do

    W_eps = integrate(this%gr%pt%x, dW_eps_dx)

    ! Finish

    return

  end function W_eps

  !****

  function tau_ss (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: tau_ss

    integer  :: j
    real(WP) :: dtau_dx(this%n)
    
    ! Evaluate the steady-state total torque, in units of G
    ! M_star**2/R_star. This expression is based on eqn. (18) of
    ! [Townsend:2018]

    !$OMP PARALLEL DO
    do j = 1, this%n
       dtau_dx(j) = this%dtau_ss_dx(j)
    end do

    tau_ss = integrate(this%gr%pt%x, dtau_dx)

    ! Finish

    return

  end function tau_ss

  !****

  function tau_tr (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: tau_tr

    integer  :: j
    real(WP) :: dtau_dx(this%n)
    
    ! Evaluate the transient total torque, in units of G
    ! M_star**2/R_star. This expression is based on eqn. (18) of
    ! [Townsend:2018]

    !$OMP PARALLEL DO
    do j = 1, this%n
       dtau_dx(j) = this%dtau_tr_dx(j)
    end do

    tau_tr = integrate(this%gr%pt%x, dtau_dx)

    ! Finish

    return

  end function tau_tr

  !****

  function omega_int (this)

    class(wave_t), intent(in) :: this
    complex(WP)               :: omega_int

    integer     :: j
    complex(WP) :: dzeta_dx(this%n)

    ! Calculate the dimensionless frequency from the zeta integral

    !OMP PARALLEL DO
    do j = 1, this%n
       dzeta_dx(j) = this%dzeta_dx(j)
    end do

    omega_int = sqrt(integrate(this%gr%pt%x, dzeta_dx))

    ! Finish

    return

  end function omega_int

  !****

  function beta (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: beta

    integer  :: j
    real(WP) :: dbeta_dx(this%n)
     
    ! Calculate the rotational splitting factor. This is based on
    ! equation (3.357) of [Aer2010]; the Ledoux constant follows from
    ! equation (3.361) [ibid] as C_nl = 1 - beta

    !$OMP PARALLEL DO
    do j = 1, this%n
       dbeta_dx(j) = this%dbeta_dx(j)
    end do

    beta = integrate(this%gr%pt%x, dbeta_dx)

    ! Finish

    return

  end function beta

  !****

  function domega_rot (this)

    class(wave_t), intent(in) :: this
    complex(WP)               :: domega_rot

    integer  :: j
    real(WP) :: dbeta_dx(this%n)
    real(WP) :: Omega_rot(this%n)

    ! Calculate the rotational splitting between modes of adjacent
    ! m. This is based on eqn. 3.355 of [Aer2010].

    !$OMP PARALLEL DO
    do j = 1, this%n

       dbeta_dx(j) = this%dbeta_dx(j)

       associate ( &
            ml => this%cx%model(), &
            pt => this%gr%pt(j))

         Omega_rot(j) = this%cx%Omega_rot(pt)

       end associate
       
    end do

    domega_rot = integrate(this%gr%pt%x, dbeta_dx*Omega_rot)

    ! Finish

    return

  end function domega_rot

  !****

  function eta (this)

    class(wave_t), intent(in) :: this
    real(WP)                  :: eta

    integer  :: j
    real(WP) :: dW_dx(this%n)
    real(WP) :: D

    ! Calculate the normalized growth rate defined (as eta') by [Stel1978]

    !$OMP PARALLEL DO
    do j = 1, this%n
       dW_dx(j) = this%dW_dx(j)
    end do

    D = integrate(this%gr%pt%x, abs(dW_dx))

    if (D /= 0._WP) then
       eta = integrate(this%gr%pt%x, dW_dx)/D
    else
       eta = 0._WP
    endif

    ! Finish

    return

  end function eta

end module gyre_wave
