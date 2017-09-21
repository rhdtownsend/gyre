! Module   : gyre_mode
! Purpose  : mode data
!
! Copyright 2013-2017 Rich Townsend
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

  use gyre_constants
  use gyre_context
  use gyre_ext
  use gyre_freq
  use gyre_grid
  use gyre_grid_util
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_state
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: mode_t
     type(c_state_t)           :: st
     type(context_t), pointer  :: cx => null()
     type(grid_t), allocatable :: gr
     type(mode_par_t)          :: md_p
     type(osc_par_t)           :: os_p
     complex(WP), allocatable  :: y_c(:,:)
     type(c_ext_t)             :: discrim
     complex(WP)               :: scl
     complex(WP)               :: omega
     complex(WP)               :: l_i
     integer                   :: n_k
     integer                   :: k_ref
     integer                   :: j
     integer                   :: l
     integer                   :: m
     integer                   :: n_pg
     integer                   :: n_p
     integer                   :: n_g
     logical                   :: pruned
   contains
     private
     procedure         :: classify_
     procedure, public :: prune
     procedure, public :: freq
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
     procedure, public :: dbeta_dx
     procedure, public :: dtau_dx_ss
     procedure, public :: dtau_dx_tr
     procedure, public :: Yt_1
     procedure, public :: Yt_2
     procedure, public :: I_0
     procedure, public :: I_1
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
     procedure, public :: H
     procedure, public :: W
     procedure, public :: W_eps
     procedure, public :: beta
     procedure, public :: omega_int
     procedure, public :: eta
  end type mode_t

  ! Interfaces

  interface mode_t
     module procedure mode_t_
  end interface mode_t

  interface reallocate
     module procedure reallocate_1_
  end interface reallocate

  $if ($MPI)
  interface bcast
     module procedure bcast_
  end interface bcast
  $endif

  ! Access specifiers

  private

  public :: mode_t
  public :: reallocate
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  function mode_t_ (st, y_c, discrim, cx, gr, md_p, os_p, j, normalize) result (md)

    type(c_state_t), intent(in)          :: st
    complex(WP), intent(in)              :: y_c(:,:)
    type(c_ext_t), intent(in)            :: discrim
    type(context_t), pointer, intent(in) :: cx
    type(grid_t), intent(in)             :: gr
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    integer, intent(in)                  :: j
    logical, intent(in), optional        :: normalize
    type(mode_t)                         :: md

    logical     :: normalize_
    real(WP)    :: x_ref
    complex(WP) :: y_1_ref
    complex(WP) :: f_phase

    $CHECK_BOUNDS(SIZE(y_c, 1),6)
    $CHECK_BOUNDS(SIZE(y_c, 2),gr%n_k)

    if (PRESENT(normalize)) then
       normalize_ = normalize
    else
       normalize_ = .TRUE.
    endif

    ! Construct the mode_t

    md%st = st
    md%cx => cx
    allocate(md%gr, SOURCE=gr)

    md%md_p = md_p
    md%os_p = os_p

    md%y_c = y_c

    md%omega = st%omega
    md%discrim = discrim

    md%l_i = cx%l_i(st)

    md%n_k = gr%n_k

    md%j = j
    md%l = md_p%l
    md%m = md_p%m
    
    md%pruned = .FALSE.

    ! Locate the reference point

    x_ref = MIN(MAX(os_p%x_ref, gr%pt(1)%x), gr%pt(gr%n_k)%x)

    md%k_ref = MINLOC(ABS(gr%pt%x - x_ref), DIM=1)

    ! Normalize so that y_1 at the reference point is purely real, and
    ! the total mode inertia E is unity

    if (normalize_) then

       md%scl = 1._WP
    
       y_1_ref = y_c(1, md%k_ref)

       f_phase = CONJG(y_1_ref)/ABS(y_1_ref)

       md%scl = 1._WP/SQRT(md%E())*f_phase

    else

       md%scl = 1._WP

    endif

    ! Classify the mode

    call md%classify_()

    ! Finish

    return

  end function mode_t_

  !****

  subroutine classify_ (this)

    class(mode_t), intent(inout) :: this

    integer  :: k
    real(WP) :: y_1(this%n_k)
    real(WP) :: y_2(this%n_k)
    integer  :: k_i
    integer  :: k_o
    real(WP) :: x_i
    integer  :: n_c
    integer  :: n_a

    ! Classify the mode based on its eigenfunctions

    if (this%l == 0) then

       ! Radial modes
       
       ! Look for the first monotonic segment in y_1 (this is to deal with
       ! noisy near-zero solutions at the origin)

       !$OMP PARALLEL DO
       do k = 1, this%n_k
          y_1(k) = REAL(this%y_i(1, k))
          y_2(k) = REAL(this%y_i(2, k))
       end do

       k_i = this%n_k

       mono_loop : do k = 2, this%n_k-1
          if ((y_1(k) >= y_1(k-1) .AND. y_1(k+1) >= y_1(k)) .OR. &
              (y_1(k) <= y_1(k-1) .AND. y_1(k+1) <= y_1(k))) then
             k_i = k
             exit mono_loop
          endif
       end do mono_loop

       ! Count winding numbers

       call count_windings_(y_1(k_i:), y_2(k_i:), n_c, n_a)

       ! Classify (the additional 1 is for the node at the center)

       this%n_p = n_a + 1
       this%n_g = n_c

       this%n_pg = this%n_p - this%n_g

    elseif (this%l == 1 .AND. .NOT. this%os_p%cowling_approx) then

       ! Dipole modes (non-Cowling)

       ! Set up the Takata Y^a_1 and Y^a_2 functions

       !$OMP PARALLEL DO
       do k = 1, this%n_k
          y_1(k) = REAL(this%Yt_1(k))
          y_2(k) = REAL(this%Yt_2(k))
       end do

       ! Find the inner turning point (this is to deal with noisy
       ! near-zero solutions at the inner boundary)

       call find_turn(this%cx%ml, this%gr, REAL(this%st%omega), this%md_p, this%os_p, k_i, x_i)

       ! Count winding numbers, taking care to avoid counting nodes at
       ! the center and surface

       if (y_1(this%n_k) == 0._WP) then
          k_o = this%n_k-1
       else
          k_o = this%n_k
       endif
       
       call count_windings_(y_1(k_i:k_o), y_2(k_i:k_o), n_c, n_a)

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
       do k = 1, this%n_k
          y_1(k) = REAL(this%y_i(1, k))
          y_2(k) = REAL(this%y_i(2, k) + this%y_i(3, k))
       end do

       ! Handle special case where the inner boundary y_1 = 0 is
       ! appled off-center -- don't count the node there

       if (this%os_p%inner_bound == 'ZERO_R') then
          k_i = 2
       else
          k_i = 1
       endif

       k_o = this%n_k

       ! Count winding numbers

       call count_windings_(y_1(k_i:k_o), y_2(k_i:k_o), n_c, n_a)

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

      integer  :: k
      real(WP) :: y_2_cross

      $CHECK_BOUNDS(SIZE(y_2),SIZE(y_1))

      if(PRESENT(x)) then
         $CHECK_BOUNDS(SIZE(x),SIZE(y_1))
      endif

      ! Count clockwise (n_c) and anticlockwise (n_a) windings in the (y_1,y_2) plane

      n_c = 0
      n_a = 0

      do k = 1,SIZE(y_1)-1

         ! Look for a node in y_1

         if (y_1(k) >= 0._WP .AND. y_1(k+1) < 0._WP) then

            ! Solve for the crossing ordinate

            y_2_cross = y_2(k) - y_1(k)*(y_2(k+1) - y_2(k))/(y_1(k+1) - y_1(k))

            if(y_2_cross >= 0._WP) then
               n_a = n_a + 1
               if(PRESENT(x)) print *,'A node:',x(k),x(k+1)
            else
               n_c = n_c + 1
               if(PRESENT(x)) print *,'C node:',x(k),x(k+1)
            endif

         elseif (y_1(k) <= 0._WP .AND. y_1(k+1) > 0._WP) then

            ! Solve for the crossing ordinate

            y_2_cross = y_2(k) - y_1(k)*(y_2(k+1) - y_2(k))/(y_1(k+1) - y_1(k))

            if (y_2_cross <= 0._WP) then
               n_a = n_a + 1
               if(PRESENT(x)) print *,'A node:',x(k),x(k+1)
            else
               n_c = n_c + 1
               if(PRESENT(x)) print *,'C node:',x(k),x(k+1)
            endif

         endif

      end do

      ! Finish

      return

    end subroutine count_windings_

  end subroutine classify_

  !****

  $REALLOCATE(type(mode_t),1)

  !****

  $if($MPI)

  subroutine bcast_ (md, root_rank, ml)

    class(mode_t), intent(inout)       :: md
    integer, intent(in)                :: root_rank
    class(model_t), intent(in), target :: ml

    ! Broadcast the mode_t

    call bcast(md%cx, root_rank, ml)

    call bcast_alloc(md%gr, root_rank)

    call bcast(md%md_p, root_rank)
    call bcast(md%os_p, root_rank)

    call bcast_alloc(md%y_c, root_rank)

    call bcast(md%discrim, root_rank)
    call bcast(md%scl, root_rank)

    call bcast(md%omega, root_rank)
    call bcast(md%l_i, root_rank)

    call bcast(md%n_k, root_rank)
    call bcast(md%k_ref, root_rank)

    call bcast(md%j, root_rank)
    call bcast(md%l, root_rank)
    call bcast(md%m, root_rank)

    call bcast(md%n_p, root_rank)
    call bcast(md%n_g, root_rank)
    call bcast(md%n_pg, root_rank)

    call bcast(md%pruned, root_rank)

    ! Finish

    return

  end subroutine bcast_

  $endif

  !****

  subroutine prune (this)

    class(mode_t), intent(inout) :: this

    ! Prune the mode_t

    if (.NOT. this%pruned) then

       deallocate(this%y_c)
       deallocate(this%gr)

       this%pruned = .TRUE.

    endif

    ! Finish

    return

  end subroutine prune

  !****

  function freq (this, freq_units, freq_frame)

    class(mode_t), intent(in)          :: this
    character(*), intent(in)           :: freq_units
    character(*), optional, intent(in) :: freq_frame
    complex(WP)                        :: freq

    ! Calculate the frequency

    associate( &
         ml => this%cx%ml, &
         pt_i => this%gr%pt(1), &
         pt_o => this%gr%pt(this%n_k) )

      if (PRESENT(freq_frame)) then
         freq = freq_from_omega(this%st%omega, ml, pt_i, pt_o, freq_units, freq_frame, this%md_p, this%os_p)
      else
         freq = freq_from_omega(this%st%omega, ml, pt_i, pt_o, freq_units, 'INERTIAL', this%md_p, this%os_p)
      endif

    end associate

    ! Finish
    
    return

  end function freq

  !****

  function y_i (this, i, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: i
    integer, intent(in)       :: k
    complex(WP)               :: y_i

    ! Evaluate y(i)

    y_i = this%scl*this%y_c(i, k)

    ! Finish

    return

  end function y_i

  !****

  function xi_r (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: xi_r

    complex(WP) :: y_1

    ! Evaluate the radial displacement perturbation, in units of
    ! R_star

    associate ( &
         pt => this%gr%pt(k), &
         l_i => this%l_i )

      y_1 = this%y_i(1, k)

      if (l_i /= 1._WP) then

         if (pt%x /= 0._WP) then
            xi_r = y_1*pt%x**(l_i-1._WP)
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

  function xi_h (this, k)
    
    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: xi_h

    complex(WP) :: y_2
    complex(WP) :: y_3
    real(WP)    :: c_1
    real(WP)    :: Omega_rot
    complex(WP) :: omega_c

    ! Evaluate the horizontal displacement perturbation, in units of
    ! R_star

    if (this%l /= 0) then

       associate ( &
            ml => this%cx%ml, &
            pt => this%gr%pt(k), &
            l_i => this%l_i )

         y_2 = this%y_i(2, k)
         y_3 = this%y_i(3, k)

         c_1 = ml%coeff(I_C_1, pt)

         Omega_rot = ml%coeff(I_OMEGA_ROT, pt)

         omega_c = this%cx%omega_c(Omega_rot, this%st)
      
         if (l_i /= 1._WP) then

            if (pt%x /= 0._WP) then
               xi_h = (y_2+y_3)*pt%x**(l_i-1._WP)/(c_1*omega_c**2)
            else
               xi_h = 0._WP
            end if

         else
            
            xi_h = (y_2+y_3)/(c_1*omega_c**2)

         endif

       end associate

    else

       xi_h = 0._WP

    end if

    ! Finish

    return
    
  end function xi_h

  !****

  function eul_phi (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: eul_phi

    complex(WP) :: y_3
    real(WP)    :: c_1
    
    ! Evaluate the Eulerian gravitational potential perturbation, in
    ! units of G M_star / R_star

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k), &
         l_i => this%l_i )

      y_3 = this%y_i(3, k)

      c_1 = ml%coeff(I_C_1, pt)

      if (l_i /= 0._WP) then

         if (pt%x /= 0._WP) then
            eul_phi = y_3*pt%x**l_i/c_1
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

  function deul_phi (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: deul_phi

    complex(WP) :: y_4
    real(WP)    :: c_1
    
    ! Evaluate the Eulerian potential gradient (gravity) perturbation,
    ! in units of G M_star / R_star**2

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k), &
         l_i => this%l_i )

      y_4 = this%y_i(4, k)

      c_1 = ml%coeff(I_C_1, pt)

      if (l_i /= 1._WP) then

         if (pt%x /= 0._WP) then
            deul_phi = y_4*pt%x**(l_i-1._WP)/c_1
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

  function lag_S (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: lag_S

    complex(WP) :: y_5

    ! Evaluate the Lagrangian specific entropy perturbation, in units
    ! of c_p

    associate ( &
         pt => this%gr%pt(k), &
         l_i => this%l_i )

      y_5 = this%y_i(5, k)

      if (pt%x /= 0._WP) then
         lag_S = y_5*pt%x**(l_i-2._WP)
      else
         lag_S = 0._WP
      endif

    end associate

    ! Finish

    return

  end function lag_S

  !****

  function lag_L (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: lag_L

    complex(WP) :: y_6

    ! Evaluate the Lagrangian radiative luminosity perturbation, in
    ! units of L_star

    associate ( &
         pt => this%gr%pt(k), &
         l_i => this%l_i )

      y_6 = this%y_i(6, k)

      if (pt%x /= 0._WP) then
         lag_L = y_6*pt%x**(l_i+1._WP)
      else
         lag_L = 0._WP
      endif

    end associate

    ! Finish

    return

  end function lag_L

  !****

  function eul_P (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: eul_P

    complex(WP) :: xi_r
    complex(WP) :: lag_P
    real(WP)    :: V_2

    ! Evaluate the Eulerian pressure perturbation, in units of P

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      xi_r = this%xi_r(k)
      lag_P = this%lag_P(k)

      V_2 = ml%coeff(I_V_2, pt)

      eul_P = lag_P + V_2*pt%x*xi_r

    end associate

    ! Finish

    return

  end function eul_P
  
  !****

  function lag_P (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: lag_P

    complex(WP) :: y_1
    complex(WP) :: y_2
    real(WP)    :: V_2

    ! Evaluate the Lagrangian pressure perturbation, in units of P

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k), &
         l_i => this%l_i )

      y_1 = this%y_i(1, k)
      y_2 = this%y_i(2, k)

      V_2 = ml%coeff(I_V_2, pt)

      if (l_i /= 0._WP) then

         if (pt%x /= 0._WP) then
            lag_P = V_2*(y_2 - y_1)*pt%x**l_i
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

  function eul_rho (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: eul_rho

    complex(WP) :: xi_r
    complex(WP) :: lag_rho
    real(WP)    :: U
    real(WP)    :: dU
    real(WP)    :: D

    ! Evaluate the Eulerian density perturbation, in units of rho

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      xi_r = this%xi_r(k)
      lag_rho = this%lag_rho(k)

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

  function lag_rho (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: lag_rho

    complex(WP) :: lag_P
    complex(WP) :: lag_S
    real(WP)    :: Gamma_1
    real(WP)    :: delta

    ! Evaluate the Lagrangian density perturbation, in units of
    ! rho. This expression implements eqn. 13.83 of [Unn1989]

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      lag_P = this%lag_P(k)
      lag_S = this%lag_S(k)

      Gamma_1 = ml%coeff(I_GAMMA_1, pt)
      delta = ml%coeff(I_DELTA, pt)

      lag_rho = lag_P/Gamma_1 - delta*lag_S

    end associate

    ! Finish

    return
    
  end function lag_rho

  !****

  function eul_T (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: eul_T

    complex(WP) :: xi_r
    complex(WP) :: lag_T
    real(WP)    :: V_2
    real(WP)    :: nabla

    ! Evaluate the Lagrangian temperature perturbation, in units of T

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      xi_r = this%xi_r(k)
      lag_T = this%lag_T(k)

      V_2 = ml%coeff(I_V_2, pt)
      nabla = ml%coeff(I_NABLA, pt)
      
      eul_T = lag_T + nabla*V_2*pt%x*xi_r

    end associate

    ! Finish

    return

  end function eul_T

  !****

  function lag_T (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: lag_T

    complex(WP) :: lag_P
    complex(WP) :: lag_S
    real(WP)    :: nabla_ad

    ! Evaluate the Lagrangian temperature perturbation, in units of
    ! T. This expression implements eqn. 13.84 of [Unn1989]

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      lag_P = this%lag_P(k)
      lag_S = this%lag_S(k)

      nabla_ad = ml%coeff(I_NABLA_AD, pt)
      
      lag_T = nabla_ad*lag_P + lag_S

    end associate

    ! Finish

    return

  end function lag_T

  !****

  function lambda (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: lambda

    real(WP) :: Omega_rot

    ! Evaluate the angular eigenvalue

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )
      
      Omega_rot = ml%coeff(I_OMEGA_ROT, pt)
    
      lambda = this%cx%lambda(Omega_rot, this%st)

    end associate

    ! Finish

    return

  end function lambda
    
  !****

  function dE_dx (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    real(WP)                  :: dE_dx

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    complex(WP) :: lambda
    real(WP)    :: U
    real(WP)    :: c_1

    ! Evaluate the differential mode inertia, in units of M_star
    ! R_star**2. This expression is based on eqn. 3.139 of [Aer2010]

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      xi_r = this%xi_r(k)
      xi_h = this%xi_h(k)

      lambda = this%lambda(k)

      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      dE_dx = (ABS(xi_r)**2 + ABS(lambda)*ABS(xi_h)**2)*U*pt%x**2/c_1

    end associate

    ! Finish

    return

  end function dE_dx

  !****

  function dW_dx (this, k)

    use gyre_evol_model

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    real(WP)                  :: dW_dx

    real(WP)    :: t_dyn
    real(WP)    :: t_kh
    complex(WP) :: lag_T
    complex(WP) :: lag_S
    real(WP)    :: c_thk

    ! Evaluate the differential work, in units of G M_star**2/R_star.
    ! This expression is based on eqn. 25.9 of [Unn1989]

    select type (ml => this%cx%ml)
    class is (evol_model_t)
       t_dyn = SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))
       t_kh = (G_GRAVITY*ml%M_star**2/ml%R_star)/ml%L_star
    class default
       t_dyn = 1._WP
       t_kh = 1._WP
    end select

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      lag_T = this%lag_T(k)
      lag_S = this%lag_S(k)
    
      c_thk = ml%coeff(I_C_THK, pt)

      dW_dx = PI*AIMAG(CONJG(lag_T)*lag_S)*c_thk*pt%x**2*t_dyn/t_kh

    end associate

    ! Finish

    return

  end function dW_dx

  !****

  function dW_eps_dx (this, k)

    use gyre_evol_model

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
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
    ! based on eqn. 25.9 of [Unn1989]
    
    select type (ml => this%cx%ml)
    class is (evol_model_t)
       t_dyn = SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star))
       t_kh = (G_GRAVITY*ml%M_star**2/ml%R_star)/ml%L_star
    class default
       t_dyn = 1._WP
       t_kh = 1._WP
    end select

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      lag_rho = this%lag_rho(k)
      lag_T = this%lag_T(k)

      c_eps = ml%coeff(I_C_EPS, pt)

      select case (this%os_p%deps_scheme)
      case ('MODEL')
         eps_rho = ml%coeff(I_EPS_RHO, pt)
         eps_T = ml%coeff(I_EPS_T, pt)
      case default
         $ABORT(Evaluating dW_dx_eps not supported for deps_scheme /= MODEL)
      end select
      
      omega_r = REAL(this%st%omega)

      dW_eps_dx = PI/omega_r*REAL(CONJG(lag_T)*c_eps*(eps_rho*lag_rho + eps_T*lag_T))*pt%x**2*t_dyn/t_kh

    end associate

    ! Finish

    return

  end function dW_eps_dx

  !****

  function dbeta_dx (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    real(WP)                  :: dbeta_dx

    complex(WP) :: xi_r
    complex(WP) :: xi_h
    real(WP)    :: U
    real(WP)    :: c_1
    complex(WP) :: lambda
    real(WP)    :: E
    
    ! Calculate the (unnormalized) rotation splitting kernel. This is
    ! based on the derivative of equation 3.357 of [Aer2010] with
    ! respect to x

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      xi_r = this%xi_r(k)
      xi_h = this%xi_h(k)

      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      ! Question: should the following be lambda or l(l+1)?

      lambda = this%lambda(k)

      E = this%E()

      dbeta_dx = 4._WP*PI*REAL((ABS(xi_r)**2 + (lambda-1._WP)*ABS(xi_h)**2 - &
                                2._WP*xi_r*CONJG(xi_h))*U*pt%x**2/c_1)/E

    end associate

    ! Finish

    return

  end function dbeta_dx

  !****

  function dtau_dx_ss (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    real(WP)                  :: dtau_dx_ss

    complex(WP) :: lag_P
    complex(WP) :: lag_rho
    real(WP)    :: V_2
    real(WP)    :: c_1
    real(WP)    :: U
    
    ! Evaluate the steady-state differential torque, in units of G
    ! M_star**2/R_star. This expression is based on eqn. 13 of
    ! [Tow2017]

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k), &
         m => this%m )

      lag_P = this%lag_P(k)

      lag_rho = this%lag_rho(k)

      V_2 = ml%coeff(I_V_2, pt)
      c_1 = ml%coeff(I_C_1, pt)
      U = ml%coeff(I_U, pt)

      dtau_dx_ss = m*pt%x**2*AIMAG(lag_rho*CONJG(lag_P))*(U/(2._WP*c_1**2*V_2))
      
    end associate

    ! Finish

    return

  end function dtau_dx_ss

  !****

  function dtau_dx_tr (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    real(WP)                  :: dtau_dx_tr

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

    ! Evaluate the steady-state differential torque, in units of G
    ! M_star**2/R_star. This expression is based on eqn. 14 of
    ! [Tow2017]

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k), &
         m => this%m )

      xi_r = this%xi_r(k)

      eul_P = this%eul_P(k)

      lag_rho = this%lag_rho(k)
      eul_rho = this%eul_rho(k)

      eul_phi = this%eul_phi(k)

      V_2 = ml%coeff(I_V_2, pt)
      c_1 = ml%coeff(I_C_1, pt)
      U = ml%coeff(I_U, pt)

      Omega_rot = ml%coeff(I_OMEGA_ROT, pt)

      omega_c = this%cx%omega_c(Omega_rot, this%st)

      dtau_dx_tr = m*pt%x**2*AIMAG((omega_c/CONJG(omega_c) - 1._WP)*( &
           lag_rho*CONJG(eul_P)/(c_1*V_2) + &
           eul_rho*CONJG(eul_phi) + &
           xi_r*CONJG(eul_rho)*pt%x/c_1))*(U/(2._WP*c_1))
           
    end associate

    ! Finish

    return
    
  end function dtau_dx_tr

  !****

  function Yt_1 (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: Yt_1

    complex(WP) :: y_1
    complex(WP) :: y_3
    complex(WP) :: y_4
    real(WP)    :: J

    ! Evaluate the Takata Y_1 function. This expression is equivalent to
    ! eqn. 69 of [Tak2006b], divided by x**(2-l)

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      y_1 = this%y_i(1, k)
      y_3 = this%y_i(3, k)
      y_4 = this%y_i(4, k)

      J = 1._WP - ml%coeff(I_U, pt)/3._WP

      Yt_1 = J*y_1 + (y_3 - y_4)/3._WP

    end associate

    ! Finish

    return

  end function Yt_1

  !****

  function Yt_2 (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: Yt_2

    complex(WP) :: y_1
    complex(WP) :: y_2

    ! Evaluate the Takata Y_2 function. This expression is equivalent to 
    ! eqn. 70 of [Tak2006b], divided by V

    y_1 = this%y_i(1, k)
    y_2 = this%y_i(2, k)

    Yt_2 = y_2 - y_1

    ! Finish

    return

  end function Yt_2

  !****

  function I_0 (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    complex(WP)               :: I_0

    complex(WP) :: y_1
    complex(WP) :: y_4
    real(WP)    :: U
    real(WP)    :: c_1
    
    ! Evaluate the I_0 integral, which should be zero for radial
    ! modes. This expression is based on eqn. 42 of [Tak2006a]

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k), &
         l_i => this%l_i )

      y_1 = this%y_i(1, k)
      y_4 = this%y_i(4, k)

      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      if (pt%x /= 0._WP) then
         I_0 = pt%x**(l_i+1._WP)*(U*y_1 + y_4)/c_1
      else
         I_0 = 0._WP
      endif

    end associate

    ! Finish

    return

  end function I_0

  !****

  function I_1 (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
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
    ! modes. This expression is based on eqn. 43 of [Tak2006a]

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k), &
         l_i => this%l_i )

      y_1 = this%y_i(1, k)
      y_2 = this%y_i(2, k)
      y_3 = this%y_i(3, k)
      y_4 = this%y_i(4, k)

      U = ml%coeff(I_U, pt)
      c_1 = ml%coeff(I_C_1, pt)

      Omega_rot = ml%coeff(I_OMEGA_ROT, pt)

      omega_c = this%cx%omega_c(Omega_rot, this%st)

      if (pt%x /= 0._WP) then
         I_1 = pt%x**(l_i+2._WP)*(c_1*omega_c**2*U*y_1 - U*y_2 + &
               (c_1*omega_c**2 - 2._WP)*y_3 + (c_1*omega_c**2 - 1._WP)*y_4)/c_1**2
      else
         I_1 = 0._WP
      endif

    end associate
      
    ! Finish

    return

  end function I_1

  !****

  function prop_type (this, k)

    class(mode_t), intent(in) :: this
    integer, intent(in)       :: k
    integer                   :: prop_type

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: Omega_rot
    real(WP) :: lambda
    real(WP) :: omega_c
    real(WP) :: g_4
    real(WP) :: g_2
    real(WP) :: g_0
    real(WP) :: gamma

    ! Set up the propagation type (0 -> evanescent, 1 -> p, -1 -> g)

    associate ( &
         ml => this%cx%ml, &
         pt => this%gr%pt(k) )

      if (ml%is_vacuum(pt)) then

         prop_type = 0

      else

         ! Calculate the discriminant gamma

         V_g = ml%coeff(I_V_2, pt)*pt%x**2/ml%coeff(I_GAMMA_1, pt)
         As = ml%coeff(I_AS, pt)
         U = ml%coeff(I_U, pt)
         c_1 = ml%coeff(I_C_1, pt)

         Omega_rot = ml%coeff(I_OMEGA_ROT, pt)

         lambda = REAL(this%lambda(k))

         omega_c = REAL(this%cx%omega_c(Omega_rot, this%st))

         g_4 = -4._WP*V_g*c_1
         g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
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

    class(mode_t), intent(in) :: this
    complex(WP)               :: lag_T_eff

    complex(WP) :: xi_r
    complex(WP) :: lag_L

    ! Calculate the effective temperature perturbation at x_ref
    ! (assumed to correspond to the photosphere), in units of
    ! T_eff. This expression is based on the standard definition of
    ! effective temperature

    associate (k => this%k_ref)

      xi_r = this%xi_r(k)
      lag_L = this%lag_L(k)

      lag_T_eff = 0.25_WP*(lag_L - 2._WP*xi_r)

    end associate

    ! Finish

    return

  end function lag_T_eff

  !****

  function lag_g_eff (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: lag_g_eff

    complex(WP) :: xi_r
    complex(WP) :: deul_phi
    real(WP)    :: c_1
    real(WP)    :: U

    ! Calculate the effective gravity perturbation at x_ref (assumed
    ! to correspond to the photosphere), in units of the gravity. This
    ! expression is based on eqn. 24 of [Dup2002]

    associate ( &
         k => this%k_ref, &
         ml => this%cx%ml, &
         pt => this%gr%pt(this%k_ref), &
         omega => this%st%omega )

      xi_r = this%xi_r(k)
      deul_phi = this%deul_phi(k)

      c_1 = ml%coeff(I_C_1, pt)
      U = ml%coeff(I_U, pt)

      lag_g_eff = (c_1/pt%x)*deul_phi + (U - (2._WP + c_1*omega**2))*xi_r/pt%x

    end associate

    ! Finish

    return

  end function lag_g_eff

  !****

  function f_T (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: f_T

    complex(WP) :: C_T

    ! Evaluate the non-adiabatic f_T parameter. This is expression is
    ! based on eqn. 5 of [Dup2003]

    associate (k => this%k_ref)

      C_T = this%lag_T_eff()/this%xi_r(k)

      f_T = ABS(C_T)

    end associate

    ! Finish

    return

  end function f_T

  !****

  function f_g (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: f_g

    complex(WP) :: C_g

    ! Evaluate the non-adiabatic f_g parameter. This is expression is
    ! based on eqn. 6 of [Dup2003]

    associate (k => this%k_ref)

      C_g = this%lag_g_eff()/this%xi_r(k)

      f_g = -ABS(C_g)

    end associate

    ! Finish

    return

  end function f_g

  !****

  function psi_T (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: psi_T

    complex(WP) :: C_T

    ! Calculate the non-adiabatic psi_T parameter, in radians. This is
    ! expression is based on eqn. 5 of [Dup2003]

    associate (k => this%k_ref)

      C_T = this%lag_T_eff()/this%xi_r(k)

      psi_T = ATAN2(AIMAG(C_T), REAL(C_T))

    end associate

    ! Finish

    return

  end function psi_T

  !****

  function psi_g (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: psi_g

    ! Calculate the non-adiabatic psi_g parameter, in radians

    psi_g = PI

    ! Finish

    return

  end function psi_g

  !****

  function E (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: E

    integer  :: k
    real(WP) :: dE_dx(this%n_k)

    ! Calculate the mode inertia, in units of M_star R_star**2

    !$OMP PARALLEL DO
    do k = 1, this%n_k
       dE_dx(k) = this%dE_dx(k)
    end do

    E = integrate(this%gr%pt%x, dE_dx)

    ! Finish

    return

  end function E

  !****

  function E_p (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: E_p

    integer  :: k
    real(WP) :: dE_dx(this%n_k)

    ! Calculate the mode inertia in acoustic-wave propagation regions,
    ! in units of M_star R_star**2

    !$OMP PARALLEL DO
    do k = 1, this%n_k
       if (this%prop_type(k) == 1) then
          dE_dx(k) = this%dE_dx(k)
       else
          dE_dx(k) = 0._WP
       endif
    end do

    E_p = integrate(this%gr%pt%x, dE_dx)

    ! Finish

    return

  end function E_p

  !****

  function E_g (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: E_g

    integer  :: k
    real(WP) :: dE_dx(this%n_k)

    ! Calculate the mode inertia in gravity-wave propagation regions,
    ! in units of M_star R_star**2

    !$OMP PARALLEL DO
    do k = 1, this%n_k
       if (this%prop_type(k) == -1) then
          dE_dx(k) = this%dE_dx(k)
       else
          dE_dx(k) = 0._WP
       endif
    end do

    E_g = integrate(this%gr%pt%x, dE_dx)

    ! Finish

    return

  end function E_g

  !****

  function E_norm (this)
 
    class(mode_t), intent(in) :: this
    real(WP)                  :: E_norm

    real(WP)    :: E
    complex(WP) :: xi_r
    complex(WP) :: xi_h
    complex(WP) :: lambda
    real(WP)    :: A2

    ! Calculate the normalized mode inertia. This expression is based
    ! on eqn. 3.140 of [Aer2010]

    associate (k => this%k_ref)

      E = this%E()

      xi_r = this%xi_r(k)
      xi_h = this%xi_h(k)

      lambda = this%lambda(k)

      select case (this%os_p%inertia_norm)
      case ('RADIAL')
         A2 = ABS(xi_r)**2
      case ('HORIZ')
         A2 = ABS(lambda)*ABS(xi_h)**2
      case ('BOTH')
         A2 = ABS(xi_r)**2 + ABS(lambda)*ABS(xi_h)**2
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

  function H (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: H

    ! Calculate the mode total energy, in units of G M_star**2 / R_star

    H = 0.5_WP*REAL(this%st%omega)**2*this%E()

    ! Finish

    return

  end function H

  !****

  function W (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: W

    integer  :: k
    real(WP) :: dW_dx(this%n_k)
    
    ! Calculate the total work, in units of G M_star**2/R_star

    !$OMP PARALLEL DO
    do k = 1, this%n_k
       dW_dx(k) = this%dW_dx(k)
    end do

    W = integrate(this%gr%pt%x, dW_dx)

    ! Finish

    return

  end function W

  !****

  function W_eps (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: W_eps
    
    integer  :: k
    real(WP) :: dW_eps_dx(this%n_k)
    
    ! Calculate the total work associated with nuclear processes, in
    ! units of G M_star**2/R_star

    !$OMP PARALLEL DO
    do k = 1, this%n_k
       dW_eps_dx(k) = this%dW_eps_dx(k)
    end do

    W_eps = integrate(this%gr%pt%x, dW_eps_dx)

    ! Finish

    return

  end function W_eps

  !****

  function beta (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: beta

    integer  :: k
    real(WP) :: dbeta_dx(this%n_k)
     
    ! Calculate the rotational splitting factor. This is based on
    ! equation 3.357 of [Aer2010]; the Ledoux constant follows from
    ! equation 3.361 [ibid] as C_nl = 1 - beta

    !$OMP PARALLEL DO
    do k = 1, this%n_k
       dbeta_dx(k) = this%dbeta_dx(k)
    end do

    beta = integrate(this%gr%pt%x, dbeta_dx)

    ! Finish

    return

  end function beta

  !****

  function omega_int (this)

    class(mode_t), intent(in) :: this
    complex(WP)               :: omega_int

    integer       :: k
    complex(WP)   :: xi_r
    complex(WP)   :: eul_phi
    complex(WP)   :: eul_rho
    complex(WP)   :: lag_rho
    complex(WP)   :: lag_P
    real(WP)      :: V_2
    real(WP)      :: As
    real(WP)      :: U
    real(WP)      :: c_1
    real(WP)      :: Gamma_1
    real(WP)      :: V_g
    real(WP)      :: x4_V
    complex(WP)   :: W_th
    complex(WP)   :: W_re
    complex(WP)   :: W_gr
    complex(WP)   :: W_xi
    complex(WP)   :: f_th(this%n_k)
    complex(WP)   :: f_re(this%n_k)
    complex(WP)   :: f_gr(this%n_k)
    complex(WP)   :: f_xi(this%n_k)

    ! Calculate the dimensionless frequency from the integral
    ! expression in eqn. (1.71) of [Dup2002]

    !$OMP PARALLEL DO PRIVATE (xi_r, eul_phi, eul_rho, lag_rho, lag_P, V_2, As, U, c_1, Gamma_1, V_g, x4_V)
    do k = 1, this%n_k

       associate ( &
            ml => this%cx%ml, &
            pt => this%gr%pt(k) )

         xi_r = this%xi_r(k)
         eul_phi = this%eul_phi(k)
         eul_rho = this%eul_rho(k)
         lag_rho = this%lag_rho(k)
         lag_P = this%lag_P(k)

         V_2 = ml%coeff(I_V_2, pt)
         As = ml%coeff(I_AS, pt)
         U = ml%coeff(I_U, pt)
         c_1 = ml%coeff(I_C_1, pt)

         Gamma_1 = ml%coeff(I_GAMMA_1, pt)

         V_g = V_2*pt%x**2/Gamma_1
         x4_V = pt%x**2/V_2

         f_th(k) = CONJG(lag_rho)*lag_P*(U*x4_V/(c_1**2))
         f_re(k) = 2._WP*REAL(lag_rho*CONJG(xi_r)*(pt%x/c_1)*(pt%x**2*U/c_1))
         f_gr(k) = CONJG(eul_rho)*eul_phi*(pt%x**2*U/c_1)
         f_xi(k) = -ABS(xi_r)**2*(pt%x/c_1)*(pt%x*U*(-V_g-As)/c_1)

       end associate

    end do

    W_th = integrate(this%gr%pt%x, f_th)
    W_re = integrate(this%gr%pt%x, f_re)
    W_gr = integrate(this%gr%pt%x, f_gr)
    W_xi = integrate(this%gr%pt%x, f_xi)

    omega_int = SQRT((W_th + W_re + W_gr + W_xi)/this%E())

    ! Finish

    return

  end function omega_int

  !****

  function eta (this)

    class(mode_t), intent(in) :: this
    real(WP)                  :: eta

    integer  :: k
    real(WP) :: dW_dx(this%n_k)
    real(WP) :: D

    ! Calculate the normalized growth rate defined (as eta') by [Stel1978]

    !$OMP PARALLEL DO
    do k = 1, this%n_k
       dW_dx(k) = this%dW_dx(k)
    end do

    D = integrate(this%gr%pt%x, ABS(dW_dx))

    if (D /= 0._WP) then
       eta = integrate(this%gr%pt%x, dW_dx)/D
    else
       eta = 0._WP
    endif

    ! Finish

    return

  end function eta

end module gyre_mode
