! Module   : gyre_grid
! Purpose  : grid construction
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

module gyre_grid

  ! Uses

  use core_kinds

  use gyre_constants
  use core_func
  use gyre_grid_par
  use core_order
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_rot
  use gyre_rot_factory
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions (used internally for root-finding)

  ! type, extends (func_t) :: geom_func_t
  !    real(WP) :: s
  !    integer  :: n
  !  contains
  !    procedure :: eval_c_ => eval_geom_func_
  ! end type geom_func_t

  type, extends (func_t) :: gamma_func_t
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     real(WP)                    :: omega
     integer                     :: s
   contains
     procedure :: eval_c_ => eval_gamma_func_
  end type gamma_func_t

  ! Access specifiers

  private

!  public :: grid_range
  public :: build_grid
!  public :: create_uniform
!  public :: create_geom
!  public :: create_log
  public :: find_turn

  ! Procedures

contains

  subroutine build_grid (ml, omega, gr_p, md_p, os_p, s, x, verbose)

    class(model_t), pointer, intent(in) :: ml
    real(WP), intent(in)                :: omega(:)
    type(grid_par_t), intent(in)        :: gr_p(:)
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    integer, allocatable, intent(out)   :: s(:)
    real(WP), allocatable, intent(out)  :: x(:)
    logical, optional, intent(in)       :: verbose

    logical :: write_info
    integer :: n
    integer :: i

    $ASSERT(SIZE(gr_p) >= 1,Empty grid_par_t)

    if (PRESENT(verbose)) then
       write_info = verbose .AND. check_log_level('INFO')
    else
       write_info = check_log_level('DEBUG')
    endif

    ! Build a grid using the supplied list of grid_par_t

    if (write_info) then

       write(OUTPUT_UNIT, 100) 'Building x grid'
100    format(A,1X,A)

    endif

    ! Create the base grid

    select case (gr_p(1)%op_type)
    case ('CREATE_CLONE')
       call ml%scaffold(s, x)
!    case ('CREATE_MIDPOINT')
!       x = [ms%x(1),0.5_WP*(ms%x(:n-1)+ms%x(2:)),ms%x(n)]
!    case ('CREATE_UNIFORM')
!       call create_uniform(gr_p(1)%n, x)
!       x = (1._WP-x)*ms%x(1) + x*ms%x(n)
!    case ('CREATE_GEOM')
!       call create_geom(gp(1)%s, gp(1)%n, x)
!       x = (1._WP-x)*ms%x(1) + x*ms%x(n)
!    case ('CREATE_LOG')
!       call create_log(gp(1)%s, gp(1)%n, x)
!       x = (1._WP-x)*ms%x(1) + x*ms%x(n)
!    case ('CREATE_FROM_FILE')
!       $ABORT(Not currently supported)
!       call create_from_file_(gp(1)%file, x)
!       x = (1._WP-x)*gp(1)%x_i + x*gp(1)%x_o
    case default
       $ABORT(Invalid op_type (the first op_type must be CREATE_*))
    end select

    if (write_info) then
       write(OUTPUT_UNIT, 110) 'x points :', SIZE(x), '[after', TRIM(gr_p(1)%op_type), 'op]'
110    format(2X,A,1X,I0,1X,A,1X,A,1X,A)
    endif

    ! Resample as necessary

    gr_p_loop : do i = 2, SIZE(gr_p)

       select case (gr_p(i)%op_type)
       case ('RESAMP_DISPERSION')
          call resample_dispersion_(ml, omega, gr_p(i), md_p, os_p, s, x)
       case ('RESAMP_THERMAL')
          call resample_thermal_(ml, omega, gr_p(i), md_p, os_p, s, x)
       case ('RESAMP_CENTER')
          call resample_center_(ml, omega, gr_p(i), md_p, os_p, s, x)
       case ('RESAMP_STRUCT')
          call resample_struct_(ml, gr_p(i), s, x)
       case ('RESAMP_UNIFORM')
          call resample_uniform_(gr_p(i), s, x)
       case default
          $ABORT(Invalid op_type (the second and subsequent op_types must be RESAMPLE_*))
       end select

       if (write_info) then
          write(OUTPUT_UNIT, 110) 'x points :', SIZE(x), '[after', TRIM(gr_p(i)%op_type), 'op]'
       endif

    end do gr_p_loop

    if (write_info) then

       write(OUTPUT_UNIT, 120) 'x range  :', MINVAL(x), '->', MAXVAL(x)
120    format(2X,A,E24.16,1X,A,1X,E24.16)

       write(OUTPUT_UNIT, *)

    endif

    ! Finish

    return

  end subroutine build_grid

  ! !****

  ! subroutine create_uniform (n, x)

  !   integer, intent(in)                :: n
  !   real(WP), allocatable, intent(out) :: x(:)

  !   integer :: i

  !   ! Create an n-point grid with uniform spacing across the [0,1]
  !   ! interval

  !   allocate(x(n))

  !   x(1) = 0._WP

  !   grid_loop : do i = 2,n-1
  !      x(i) = (i-1._WP)/(n-1._WP)
  !   end do grid_loop

  !   x(n) = 1._WP

  !   ! Finish

  !   return

  ! end subroutine create_uniform

  !****

!   subroutine create_geom (s, n, x)

!     real(WP), intent(in)               :: s
!     integer, intent(in)                :: n
!     real(WP), allocatable, intent(out) :: x(:)

!     integer           :: m
!     real(WP)          :: g_a
!     real(WP)          :: g_b
!     real(WP)          :: g
!     type(geom_func_t) :: gf
!     real(WP)          :: dx
!     integer           :: k

!     ! Create an n-point grid with geometric spacing in each half of
!     ! the [0,1] interval. The parameter s controls the ratio between
!     ! the boundary cell size and the average cell size 1/(n-1)

!     allocate(x(n))

!     if(MOD(n, 2) == 0) then

!        ! Even number of grid points / odd number of cells

!        ! Solve for the growth factor g. The upper bound is derived
!        ! by applying a Taylor expansion to the equation for g

!        m = n/2-1

!        gf%n = n
!        gf%s = s

!        g_a = EPSILON(0._WP)
!        g_b = (s*(n-1)-2*m-1)/m

!        g = gf%root(g_a, g_b, 0._WP)

!        ! Set up the inner part of the grid

!        x(1) = 0._WP
!        dx = 1._WP/(s*(n-1))
       
!        even_grid_loop : do k = 1,m
!           x(k+1) = x(k) + dx
!           dx = (1._WP+g)*dx
!        end do even_grid_loop

!        ! Reflect to get the outer part of the grid

!        x(m+2:) = 1._WP - x(m+1:1:-1)

!     else

!        ! Odd number of grid points / even number of cells

!        ! Solve for the growth factor g. The upper bound is derived
!        ! by applying a Taylor expansion to the equation for g

!        m = (n-1)/2

!        gf%n = n
!        gf%s = s

!        g_a = EPSILON(0._WP)
!        g_b = (s*(n-1)-2*m)/(m*(m-1))

!        g = gf%root(g_a, g_b, 0._WP)

!        ! Set up the inner part of the grid

!        x(1) = 0._WP
!        dx = 1._WP/(s*(n-1))
       
!        odd_grid_loop : do k = 1,m-1
!           x(k+1) = x(k) + dx
!           dx = (1._WP+g)*dx
!        end do odd_grid_loop

!        ! Reflect to get the outer part of the grid

!        x(m+1:) = 0.5_WP

!        x(m+2:) = 1._WP - x(m:1:-1)

!     end if

!     ! Finish

!     return

!   end subroutine create_geom

! !****

!   function eval_geom_func_ (this, z) result (f_z)

!     class(geom_func_t), intent(inout) :: this
!     complex(WP), intent(in)           :: z
!     complex(WP)                       :: f_z

!     real(WP) :: g
!     integer  :: m

!     ! Calcuate the discriminant for the geom grid growth factor

!     g = REAL(z)

!     if(MOD(this%n, 2) == 0) then

!        m = this%n/2-1

!        if(1._WP+g > HUGE(0._WP)**(1._WP/m)) then
!           f_z = - (2._WP + g)
!        else
!           f_z = (2._WP + this%s*(this%n-1)*g)/(1._WP + g)**m - (2._WP + g)
!        endif

!     else

!        m = (this%n-1)/2

!        if(1._WP+g > HUGE(0._WP)**(1._WP/m)) then
!           f_z = -2._WP
!        else
!           f_z = (2._WP + this%s*(this%n-1)*g)/(1._WP + g)**m - 2._WP
!        endif

!     endif

!     ! Finish

!     return

!   end function eval_geom_func_

! !****

!   subroutine create_log (s, n, x)

!     real(WP), intent(in)               :: s
!     integer, intent(in)                :: n
!     real(WP), allocatable, intent(out) :: x(:)

!     real(WP) :: dx_1
!     integer  :: k
!     real(WP) :: w
!     real(WP) :: t

!     ! Create an n-point grid with logarithmic spacing in each half of
!     ! the [0,1] interval. The parameter s controls the ratio between
!     ! the mean cell size and the boundary cell size. omega is not used

!     allocate(x(n))

!     dx_1 = 1._WP/(s*(n-1))

!     if (MOD(n, 2) == 0) then

!        ! Even number of grid points / odd number of cells

!        ! Set up the inner part of the grid

!        x(1) = 0._WP

!        even_grid_loop : do k = 2,n/2

!           w = (k-1.5_WP)/(n/2-1.5_WP)
!           t = (1._WP-w)*LOG(0.5_WP) + w*LOG(dx_1)

!           x(n/2-k+2) = EXP(t)

!        enddo even_grid_loop
       
!        ! Reflect to get the outer part of the grid

!        x(n/2+1:) = 1._WP - x(n/2:1:-1)

!     else

!        ! Odd number of grid points / even number of cells

!        ! Set up the inner part of the grid

!        x(1) = 0._WP

!        odd_grid_loop : do k = 2,(n-1)/2

!           w = (k-1._WP)/((n-1)/2-1._WP)
!           t = (1._WP-w)*LOG(0.5_WP) + w*LOG(dx_1)

!           x((n-1)/2-k+2) = EXP(t)

!        end do odd_grid_loop

!        x((n+1)/2) = 0.5_WP

!        ! Reflect to get the outer part of the grid

!        x((n+1)/2+1:) = 1._WP - x((n-1)/2:1:-1)

!     end if

!     ! Finish

!     return

!   end subroutine create_log

! !****

!   subroutine create_from_file_ (file, x)

!     character(LEN=*), intent(in)       :: file
!     real(WP), allocatable, intent(out) :: x(:)

!     integer :: unit
!     integer :: n
!     integer :: k

!     ! Create a grid by reading from the file

!     ! Count lines

!     open(NEWUNIT=unit, FILE=file, STATUS='OLD')

!     n = 0

!     count_loop : do
!        read(unit, *, END=100)
!        n = n + 1
!     end do count_loop

! 100 continue

!     ! Read data

!     rewind(unit)

!     allocate(x(n))

!     read_loop : do k = 1,n
!        read(unit, *) x(k)
!     end do read_loop

!     close(unit)

!     $ASSERT(x(1)==0._WP,First point not at zero)
!     $ASSERT(x(n)==1._WP,Last point not at one)

!     ! Finish

!     return

!   end subroutine create_from_file_

!****

  subroutine resample_ (dn, s, x)

    integer, intent(in)                  :: dn(:)
    integer, allocatable, intent(inout)  :: s(:)
    real(WP), allocatable, intent(inout) :: x(:)

    integer               :: n_k
    integer, allocatable  :: s_new(:)
    real(WP), allocatable :: x_new(:)
    integer               :: k_new
    integer               :: k
    integer               :: i
    
    $CHECK_BOUNDS(SIZE(s),SIZE(dn)+1)
    $CHECK_BOUNDS(SIZE(x),SIZE(dn)+1)

    ! Resample the grid by inserting dn additional points on a
    ! cell-by-cell basis
    
    n_k = SIZE(s)

    allocate(s_new(n_k + SUM(dn)))
    allocate(x_new(n_k + SUM(dn)))

    k_new = 1

    do k = 1, n_k-1

       do i = 1, dn(k)+1

          s_new(k_new) = s(k)
          x_new(k_new) = x(k) + (i-1)*(x(k+1)-x(k))/(dn(k)+1)

          k_new = k_new + 1

       end do

    end do
    
    s_new(k_new) = s(n_k)
    x_new(k_new) = x(n_k)

    call MOVE_ALLOC(s_new, s)
    call MOVE_ALLOC(x_new, x)

    ! Finish

    return

  end subroutine resample_

  !****

  subroutine resample_dispersion_ (ml, omega, gr_p, md_p, os_p, s, x)

    class(model_t), pointer, intent(in)  :: ml
    real(WP), intent(in)                 :: omega(:)
    type(grid_par_t), intent(in)         :: gr_p
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    integer, allocatable, intent(inout)  :: s(:)
    real(WP), allocatable, intent(inout) :: x(:)

    class(r_rot_t), allocatable :: rt
    integer                     :: n_k
    real(WP), allocatable       :: beta_r_max(:)
    real(WP), allocatable       :: beta_i_max(:)
    integer                     :: k
    real(WP)                    :: V_g
    real(WP)                    :: As
    real(WP)                    :: U
    real(WP)                    :: c_1
    integer                     :: j
    real(WP)                    :: omega_c
    real(WP)                    :: lambda
    real(WP)                    :: l_i
    real(WP)                    :: g_4
    real(WP)                    :: g_2
    real(WP)                    :: g_0
    real(WP)                    :: gamma
    integer, allocatable        :: dn(:)
    real(WP)                    :: dphi_osc
    real(WP)                    :: dphi_exp

    $CHECK_BOUNDS(SIZE(x),SIZE(s))

    ! Resample the grid (s,x) by adding points to each cell, such that
    ! there are at least alpha_osc points per oscillatory wavelength
    ! and alpha_exp points per exponential wavelength.
    !
    ! Wavelengths are calculated based on a local dispersion analysis
    ! of the adibatic/Cowling wave equation, for inertial frequencies
    ! specified by omega

    allocate(rt, SOURCE=r_rot_t(ml, md_p, os_p))

    ! At each point of x, determine the maximum absolute value of the
    ! real and imaginary parts of the local radial wavenumber beta

    n_k = SIZE(s)

    allocate(beta_r_max(n_k))
    allocate(beta_i_max(n_k))

    beta_r_max(1) = 0._WP
    beta_i_max(1) = 0._WP

    wavenumber_loop : do k = 2, n_k-1

       V_g = ml%V_2(s(k), x(k))*x(k)**2/ml%Gamma_1(s(k), x(k))
       As = ml%As(s(k), x(k))
       U = ml%U(s(k), x(k))
       c_1 = ml%c_1(s(k), x(k))

       beta_r_max(k) = 0._WP
       beta_i_max(k) = 0._WP

       omega_loop : do j = 1, SIZE(omega)

          omega_c = rt%omega_c(s(k), x(k), omega(j))

          lambda = rt%lambda(s(k), x(k), omega(j))
          l_i = rt%l_i(omega(j))
            
          ! Calculate the propagation discriminant gamma

          g_4 = -4._WP*V_g*c_1
          g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
          g_0 = -4._WP*lambda*As/c_1

          gamma = (g_4*omega_c**4 + g_2*omega_c**2 + g_0)/omega_c**2

          ! Update the wavenumber maxima

          if (gamma < 0._WP) then
               
             ! Propagation zone

             beta_r_max(k) = MAX(beta_r_max(k), ABS(0.5_WP*SQRT(-gamma))/x(k))
             beta_i_max(k) = MAX(beta_i_max(k), ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i))/x(k))

          else

             ! Evanescent zone

             beta_i_max(k) = MAX(beta_i_max(k), &
                  ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i - SQRT(gamma)))/x(k), &
                  ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_i + SQRT(gamma)))/x(k))
             
          end if

       end do omega_loop

    end do wavenumber_loop

    beta_r_max(n_k) = 0._WP
    beta_i_max(n_k) = 0._WP

    ! Determine how many points to add to each cell

    allocate(dn(n_k-1))

    cell_loop : do k = 1, n_k-1

       if (s(k) == s(k+1)) then

          ! Calculate the oscillatory and exponential phase change across
          ! the cell

          dphi_osc = MAX(beta_r_max(k), beta_r_max(k+1))*(x(k+1) - x(k))
          dphi_exp = MAX(beta_i_max(k), beta_i_max(k+1))*(x(k+1) - x(k))

          ! Set up dn

          dn(k) = MAX(FLOOR((gr_p%alpha_osc*dphi_osc)/TWOPI), FLOOR((gr_p%alpha_exp*dphi_exp)/TWOPI))

       else

          ! Don't add points between segments

          dn(k) = 0

       endif

    end do cell_loop

    ! Perform the resampling

    call resample_(dn, s, x)

    ! Finish

    return

  end subroutine resample_dispersion_

  !****

  subroutine resample_thermal_ (ml, omega, gr_p, md_p, os_p, s, x)

    class(model_t), pointer, intent(in)  :: ml
    real(WP), intent(in)                 :: omega(:)
    type(grid_par_t), intent(in)         :: gr_p
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    integer, allocatable, intent(inout)  :: s(:)
    real(WP), allocatable, intent(inout) :: x(:)

    class(r_rot_t), allocatable :: rt
    integer                     :: n_k
    real(WP), allocatable       :: beta_t_max(:)
    integer                     :: k
    real(WP)                    :: V
    real(WP)                    :: nabla
    real(WP)                    :: c_rad
    real(WP)                    :: c_thm
    integer                     :: j
    real(WP)                    :: omega_c
    integer, allocatable        :: dn(:)
    real(WP)                    :: dphi_thm

    $CHECK_BOUNDS(SIZE(x),SIZE(s))

    ! Resample the grid (s,x) by adding points to each cell, such that
    ! there are at least alpha_thm points per thermal length c_thm

    allocate(rt, SOURCE=r_rot_t(ml, md_p, os_p))

    ! At each point of x, determine the maximum absolute value of the
    ! local thermal wavenumber beta_t

    n_k = SIZE(s)

    allocate(beta_t_max(n_k))

    beta_t_max(1) = 0._WP

    wavenumber_loop : do k = 2, n_k-1

       V = ml%V_2(s(k), x(k))*x(k)**2
       nabla = ml%nabla(s(k), x(k))

       c_rad = ml%c_rad(s(k), x(k))
       c_thm = ml%c_thm(s(k), x(k))

       beta_t_max(k) = 0._WP

       omega_loop : do j = 1, SIZE(omega)

          omega_c = rt%omega_c(s(k), x(k), omega(j))

          beta_t_max(k) = MAX(beta_t_max(k), SQRT(ABS(V*nabla*omega_c*c_thm/c_rad))/x(k))

       end do omega_loop

    end do wavenumber_loop

    beta_t_max(n_k) = 0._WP

    ! Determine how many points to add to each cell

    allocate(dn(n_k-1))

    cell_loop : do k = 1, n_k-1

       if (s(k) == s(k+1)) then

          ! Calculate the thermal phase change across the cell

          dphi_thm = MAX(beta_t_max(k), beta_t_max(k+1))*(x(k+1) - x(k))

          ! Set up dn

          dn(k) = FLOOR((gr_p%alpha_thm*dphi_thm)/TWOPI)

       else

          ! Don't add points between segments

          dn(k) = 0

       endif

    end do cell_loop

    ! Perform the resampling

    call resample_(dn, s, x)

    ! Finish

    return

  end subroutine resample_thermal_

  !****

  subroutine resample_center_ (ml, omega, gr_p, md_p, os_p, s, x)

    class(model_t), pointer,  intent(in) :: ml
    real(WP), intent(in)                 :: omega(:)
    type(grid_par_t), intent(in)         :: gr_p
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    integer, allocatable, intent(inout)  :: s(:)
    real(WP), allocatable, intent(inout) :: x(:)
    
    integer              :: k_turn
    real(WP)             :: x_turn
    integer              :: j
    integer              :: k_turn_j
    real(WP)             :: x_turn_j
    integer              :: n_add
    integer              :: n_k
    integer, allocatable :: dn(:)
    integer              :: k

    $ASSERT(ALLOCATED(x),No input grid)

    ! Resample the grid (s,x) by adding points to central cells, such
    ! that there are n_floor points covering the evanescent region at
    ! the center. Evanescence is determined based on a local
    ! dispersion analysis of the adibatic/Cowling wave equation, for
    ! frequencies in the interval [omega_min,omega_max]

    ! First, locate the innermost turning point

    k_turn = 0
    x_turn = HUGE(0._WP)

    omega_loop : do j = 1, SIZE(omega)
       call find_turn(ml, omega(j), s, x, md_p, os_p, k_turn_j, x_turn_j)
       if (x_turn_j < x_turn) then
          k_turn = k_turn_j
          x_turn = x_turn_j
       endif
    end do omega_loop

    ! Determine how many points need to be added

    n_add = MAX(gr_p%n-k_turn, 0)

    ! Determine how many points to add to each cell

    n_k = SIZE(s)

    allocate(dn(n_k-1))

    dn = 0

    if (k_turn >= 1 .AND. k_turn < n_k) then

       do k = 1, k_turn
          if (s(k) == s(k+1)) then
             dn(k) = CEILING(n_add*(x(k_turn+1)/x_turn)/k_turn)
          endif
       end do

    end if
    
    ! Perform the resampling

    call resample_(dn, s, x)

    ! Finish

    return

  end subroutine resample_center_

  !****

  subroutine resample_struct_ (ml, gr_p, s, x)

    class(model_t), pointer, intent(in)  :: ml
    type(grid_par_t), intent(in)         :: gr_p
    integer, allocatable, intent(inout)  :: s(:)
    real(WP), allocatable, intent(inout) :: x(:)

    integer :: dn(SIZE(x)-1)
    integer :: k

    ! Resample the grid (s,x) by adding points to each cell, such that
    ! there are at least alpha_str points per dex change in the
    ! structure variables (V, As, Gamma_1, c_1, & U)

    ! Calculate the number of points to add to each cell

    dn = 0

    cell_loop : do k = 1, SIZE(x)-1

       if (s(k) == s(k+1)) then

          dn(k) = dn(k) + FLOOR(gr_p%alpha_str*dlog_(ml%V_2(s(k), x(k)), ml%V_2(s(k+1), x(k+1)))) + &
                          FLOOR(gr_p%alpha_str*dlog_(ml%As(s(k), x(k)), ml%As(s(k+1), x(k+1)))) + &
                          FLOOR(gr_p%alpha_str*dlog_(ml%Gamma_1(s(k), x(k)), ml%Gamma_1(s(k+1), x(k+1)))) + &
                          FLOOR(gr_p%alpha_str*dlog_(ml%c_1(s(k), x(k)), ml%c_1(s(k+1), x(k+1)))) + &
                          FLOOR(gr_p%alpha_str*dlog_(ml%U(s(k), x(k)), ml%U(s(k+1), x(k+1))))

       else

          dn(k) = 0

       endif

    end do cell_loop

    ! Perform the resampling

    call resample_(dn, s, x)

    ! Finish

    return

  contains

    function dlog_ (y_a, y_b) result (dlog)

      real(WP), intent(in) :: y_a
      real(WP), intent(in) :: y_b
      real(WP)             :: dlog

      ! Calculate the logarithmic change between y_a and y_b

      if ((y_a > 0._WP .AND. y_b > 0._WP) .OR. &
          (y_a < 0._WP .AND. y_b < 0._WP)) then
         dlog = ABS(LOG10(y_b/y_a))
      else
         dlog = 0._WP
      endif

      ! Finish

      return

    end function dlog_

  end subroutine resample_struct_

  !****

  subroutine resample_uniform_ (gr_p, s, x)

    type(grid_par_t), intent(in)         :: gr_p
    integer, allocatable, intent(inout)  :: s(:)
    real(WP), allocatable, intent(inout) :: x(:)

    integer              :: n_k
    integer, allocatable :: dn(:)
    integer              :: k

    $ASSERT(ALLOCATED(x),No input grid)

    ! Resample x by adding n points to each cell

    n_k = SIZE(s)

    allocate(dn(n_k-1))

    cell_loop : do k = 1, n_k-1
       if (s(k) == s(k+1)) then
          dn = gr_p%n
       else
          dn = 0
       endif
    end do cell_loop

    ! Perform the resampling

    call resample_(dn, s, x)

    ! Finish

    return

  end subroutine resample_uniform_

  !****

  subroutine find_turn (ml, omega, s, x, md_p, os_p, k_turn, x_turn)

    class(model_t), pointer, intent(in) :: ml
    real(WP), intent(in)                :: omega
    integer, intent(in)                 :: s(:)
    real(WP), intent(in)                :: x(:)
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    integer, intent(out)                :: k_turn
    real(WP), intent(out)               :: x_turn

    type(gamma_func_t) :: gf
    real(WP)           :: gf_a
    real(WP)           :: gf_b
    integer            :: k

    ! Find the cell index and location of the inner turning point at
    ! frequency omega

    x_turn = HUGE(0._WP)

    gf%ml => ml
    allocate(gf%rt, SOURCE=r_rot_t(ml, md_p, os_p))
    gf%omega = omega

    gf%s = s(1)
    gf_b = gf%eval(x(1))

    turn_loop : do k = 1, SIZE(x)-2

       gf_a = gf_b

       gf%s = s(k+1)
       gf_b = gf%eval(x(k+1))
       
       if (gf_a > 0._WP .AND. gf_b <= 0._WP) then

          k_turn = k

          if (s(k) == s(k+1)) then

             if (ABS(gf_a) < EPSILON(0._WP)*ABS(gf_b)) then
                x_turn = x(k_turn)
             elseif (ABS(gf_b) < EPSILON(0._WP)*ABS(gf_a)) then
                x_turn = x(k_turn+1)
             else
                x_turn = gf%root(x(k), x(k+1), 0._WP)
             endif

          else

             x_turn = x(k_turn)

          end if

          exit turn_loop

       endif

    end do turn_loop

    ! Finish

    return

  end subroutine find_turn

  !****

  function eval_gamma_func_ (this, z) result (gamma)

    class(gamma_func_t), intent(inout) :: this
    complex(WP), intent(in)            :: z
    complex(WP)                        :: gamma

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: g_4
    real(WP) :: g_2
    real(WP) :: g_0

    ! Calculate the propagation discriminant

    associate (s => this%s, &
               x => REAL(z))

      V_g = this%ml%V_2(s, x)*x**2/this%ml%Gamma_1(s, x)
      As = this%ml%As(s, x)
      U = this%ml%U(s, x)
      c_1 = this%ml%c_1(s, x)

      lambda = this%rt%lambda(s, x, this%omega)

      g_4 = -4._WP*V_g*c_1
      g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
      g_0 = -4._WP*lambda*As/c_1

      gamma = (g_4*this%omega**4 + g_2*this%omega**2 + g_0)/this%omega**2

    end associate

    ! Finish

    return

  end function eval_gamma_func_

end module gyre_grid
