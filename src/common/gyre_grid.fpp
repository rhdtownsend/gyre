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
  use core_order

  use gyre_grid_par
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

  type, extends (func_t) :: geom_func_t
     real(WP) :: s
     integer  :: n
   contains
     procedure :: eval_c_ => eval_geom_func_
  end type geom_func_t

  type, extends (func_t) :: gamma_func_t
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     real(WP)                    :: omega
   contains
     procedure :: eval_c_ => eval_gamma_func_
  end type gamma_func_t

  ! Access specifiers

  private

  public :: build_grid
  public :: create_uniform
  public :: create_geom
  public :: create_log
  public :: find_x_turn

  ! Procedures

contains

  subroutine build_grid (gp, ml, mp, op, omega, x_in, x, verbose)

    type(grid_par_t), intent(in)        :: gp(:)
    class(model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)        :: mp
    type(osc_par_t), intent(in)         :: op
    real(WP), intent(in)                :: omega(:)
    real(WP), allocatable, intent(in)   :: x_in(:)
    real(WP), allocatable, intent(out)  :: x(:)
    logical, optional, intent(in)       :: verbose

    logical :: write_info
    integer :: n_in
    integer :: i

    $ASSERT(SIZE(gp) >= 1,Empty grid_par_t)

    if(PRESENT(verbose)) then
       write_info = verbose .AND. check_log_level('INFO')
    else
       write_info = check_log_level('DEBUG')
    endif

    ! Build a grid using the supplied list of grid_par_t

    if(write_info) then

       write(OUTPUT_UNIT, 100) 'Building x grid'
100    format(A,1X,A)

    endif

    select case (gp(1)%op_type)
    case ('CREATE_CLONE')
       $ASSERT(ALLOCATED(x_in),No grid to clone)
       x = x_in
    case ('CREATE_MIDPOINT')
       $ASSERT(ALLOCATED(x_in),No grid to clone)
       n_in = SIZE(x_in)
       x = [x_in(1),0.5_WP*(x_in(:n_in-1)+x_in(2:)),x_in(n_in)]
    case ('CREATE_UNIFORM')
       call create_uniform(gp(1)%n, x)
    case ('CREATE_GEOM')
       call create_geom(gp(1)%s, gp(1)%n, x)
    case ('CREATE_LOG')
       call create_log(gp(1)%s, gp(1)%n, x)
    case ('CREATE_FROM_FILE')
       call create_from_file_(gp(1)%file, x)
    case default
       $ABORT(Invalid op_type (the first op_type must be CREATE_*))
    end select

    if(write_info) then
          write(OUTPUT_UNIT, 110) 'x points :', SIZE(x), '[after', TRIM(gp(1)%op_type), 'op]'
110    format(2X,A,1X,I0,1X,A,1X,A,1X,A)
    endif

    gp_loop : do i = 2,SIZE(gp)

       select case (gp(i)%op_type)
       case ('RESAMP_DISPERSION')
          call resample_dispersion_(ml, mp, op, omega, gp(i)%alpha_osc, gp(i)%alpha_exp, x)
       case ('RESAMP_THERMAL')
          call resample_thermal_(ml, mp, op, omega, gp(i)%alpha_thm, x)
       case ('RESAMP_CENTER')
          call resample_center_(ml, mp, op, omega, gp(i)%n, x)
       case ('RESAMP_STRUCT')
          call resample_struct_(ml, gp(i)%alpha_str, x)
       case ('RESAMP_UNIFORM')
          call resample_uniform_(gp(i)%n, x)
       case default
          $ABORT(Invalid op_type (the second and subsequent op_types must be RESAMPLE_*))
       end select

       if(write_info) then
          write(OUTPUT_UNIT, 110) 'x points :', SIZE(x), '[after', TRIM(gp(i)%op_type), 'op]'
       endif

    end do gp_loop

    if(write_info) then

       write(OUTPUT_UNIT, 120) 'x range  :', MINVAL(x), '->', MAXVAL(x)
120    format(2X,A,E24.16,1X,A,1X,E24.16)

       write(OUTPUT_UNIT, *)

    endif

    ! Finish

    return

  end subroutine build_grid

!****

  subroutine create_uniform (n, x)

    integer, intent(in)                :: n
    real(WP), allocatable, intent(out) :: x(:)

    integer :: i

    ! Create an n-point grid with uniform spacing across the [0,1]
    ! interval

    allocate(x(n))

    x(1) = 0._WP

    grid_loop : do i = 2,n-1
       x(i) = (i-1._WP)/(n-1._WP)
    end do grid_loop

    x(n) = 1._WP

    ! Finish

    return

  end subroutine create_uniform

!****

  subroutine create_geom (s, n, x)

    real(WP), intent(in)               :: s
    integer, intent(in)                :: n
    real(WP), allocatable, intent(out) :: x(:)

    integer           :: m
    real(WP)          :: g_a
    real(WP)          :: g_b
    real(WP)          :: g
    type(geom_func_t) :: gf
    real(WP)          :: dx
    integer           :: k

    ! Create an n-point grid with geometric spacing in each half of the
    ! [0,1] interval. The parameter s controls the ratio between the
    ! boundary cell size and the average cell size 1/(n-1)

    allocate(x(n))

    if(MOD(n, 2) == 0) then

       ! Even number of grid points / odd number of cells

       ! Solve for the growth factor g. The upper bound is derived
       ! by applying a Taylor expansion to the equation for g

       m = n/2-1

       gf%n = n
       gf%s = s

       g_a = EPSILON(0._WP)
       g_b = (s*(n-1)-2*m-1)/m

       g = gf%root(g_a, g_b, 0._WP)

       ! Set up the inner part of the grid

       x(1) = 0._WP
       dx = 1._WP/(s*(n-1))
       
       even_grid_loop : do k = 1,m
          x(k+1) = x(k) + dx
          dx = (1._WP+g)*dx
       end do even_grid_loop

       ! Reflect to get the outer part of the grid

       x(m+2:) = 1._WP - x(m+1:1:-1)

    else

       ! Odd number of grid points / even number of cells

       ! Solve for the growth factor g. The upper bound is derived
       ! by applying a Taylor expansion to the equation for g

       m = (n-1)/2

       gf%n = n
       gf%s = s

       g_a = EPSILON(0._WP)
       g_b = (s*(n-1)-2*m)/(m*(m-1))

       g = gf%root(g_a, g_b, 0._WP)

       ! Set up the inner part of the grid

       x(1) = 0._WP
       dx = 1._WP/(s*(n-1))
       
       odd_grid_loop : do k = 1,m-1
          x(k+1) = x(k) + dx
          dx = (1._WP+g)*dx
       end do odd_grid_loop

       ! Reflect to get the outer part of the grid

       x(m+1:) = 0.5_WP

       x(m+2:) = 1._WP - x(m:1:-1)

    end if

    ! Finish

    return

  end subroutine create_geom

!****

  function eval_geom_func_ (this, z) result (f_z)

    class(geom_func_t), intent(inout) :: this
    complex(WP), intent(in)           :: z
    complex(WP)                       :: f_z

    real(WP) :: g
    integer  :: m

    ! Calcuate the discriminant for the geom grid growth factor

    g = REAL(z)

    if(MOD(this%n, 2) == 0) then

       m = this%n/2-1

       if(1._WP+g > HUGE(0._WP)**(1._WP/m)) then
          f_z = - (2._WP + g)
       else
          f_z = (2._WP + this%s*(this%n-1)*g)/(1._WP + g)**m - (2._WP + g)
       endif

    else

       m = (this%n-1)/2

       if(1._WP+g > HUGE(0._WP)**(1._WP/m)) then
          f_z = -2._WP
       else
          f_z = (2._WP + this%s*(this%n-1)*g)/(1._WP + g)**m - 2._WP
       endif

    endif

    ! Finish

    return

  end function eval_geom_func_

!****

  subroutine create_log (s, n, x)

    real(WP), intent(in)               :: s
    integer, intent(in)                :: n
    real(WP), allocatable, intent(out) :: x(:)

    real(WP) :: dx_1
    integer  :: k
    real(WP) :: w
    real(WP) :: t

    ! Create an n-point grid with logarithmic spacing in each half of
    ! the [0,1] interval. The parameter s controls the ratio between
    ! the mean cell size and the boundary cell size. omega is not used

    allocate(x(n))

    dx_1 = 1._WP/(s*(n-1))

    if(MOD(n, 2) == 0) then

       ! Even number of grid points / odd number of cells

       ! Set up the inner part of the grid

       x(1) = 0._WP

       even_grid_loop : do k = 2,n/2

          w = (k-1.5_WP)/(n/2-1.5_WP)
          t = LOG(0.5_WP)*(1._WP-w) + LOG(dx_1)*w

          x(n/2-k+2) = EXP(t)

       enddo even_grid_loop
       
       ! Reflect to get the outer part of the grid

       x(n/2+1:) = 1._WP - x(n/2:1:-1)

    else

       ! Odd number of grid points / even number of cells

       ! Set up the inner part of the grid

       x(1) = 0._WP

       odd_grid_loop : do k = 2,(n-1)/2

          w = (k-1._WP)/((n-1)/2-1._WP)
          t = LOG(0.5_WP)*(1._WP-w) + LOG(dx_1)*w

          x((n-1)/2-k+2) = EXP(t)

       end do odd_grid_loop

       x((n+1)/2) = 0.5_WP

       ! Reflect to get the outer part of the grid

       x((n+1)/2+1:) = 1._WP - x((n-1)/2:1:-1)

    end if

    ! Finish

    return

  end subroutine create_log

!****

  subroutine create_from_file_ (file, x)

    character(LEN=*), intent(in)       :: file
    real(WP), allocatable, intent(out) :: x(:)

    integer :: unit
    integer :: n
    integer :: k

    ! Create a grid by reading from the file

    ! Count lines

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    n = 0

    count_loop : do
       read(unit, *, END=100)
       n = n + 1
    end do count_loop

100 continue

    ! Read data

    rewind(unit)

    allocate(x(n))

    read_loop : do k = 1,n
       read(unit, *) x(k)
    end do read_loop

    close(unit)

    ! Finish

    return

  end subroutine create_from_file_

!****

  subroutine resample_ (x, dn)

    real(WP), allocatable, intent(inout) :: x(:)
    integer, intent(in)                  :: dn(:)

    integer               :: n
    real(WP), allocatable :: x_new(:)
    integer               :: i
    integer               :: j
    integer               :: k
    
    $CHECK_BOUNDS(SIZE(dn),SIZE(x)-1)

    ! Resample x by inserting dn additional points across the cells
    
    n = SIZE(x)

    allocate(x_new(SUM(dn) + n))

    k = 1

    do i = 1,n-1
       do j = 1,dn(i)+1
          x_new(k) = x(i) + (j-1)*(x(i+1)-x(i))/(dn(i)+1)
          k = k + 1
       end do
    end do
    
    x_new(k) = x(n)

    call MOVE_ALLOC(x_new, x)

    ! Finish

    return

  end subroutine resample_

!****

  subroutine resample_dispersion_ (ml, mp, op, omega, alpha_osc, alpha_exp, x)

    class(model_t), pointer, intent(in)  :: ml
    type(mode_par_t), intent(in)         :: mp
    type(osc_par_t), intent(in)          :: op
    real(WP), intent(in)                 :: omega(:)
    real(WP), intent(in)                 :: alpha_osc
    real(WP), intent(in)                 :: alpha_exp
    real(WP), allocatable, intent(inout) :: x(:)

    class(r_rot_t), allocatable :: rt
    integer                     :: n_x
    real(WP), allocatable       :: k_r_max(:)
    real(WP), allocatable       :: k_i_max(:)
    integer                     :: i
    integer                     :: j
    real(WP)                    :: omega_c
    real(WP)                    :: lambda
    real(WP)                    :: l_0
    real(WP)                    :: g_4
    real(WP)                    :: g_2
    real(WP)                    :: g_0
    real(WP)                    :: gamma
    integer, allocatable        :: dn(:)
    real(WP)                    :: dphi_osc
    real(WP)                    :: dphi_exp

    $ASSERT(ALLOCATED(x),No input grid)

    ! Resample x by adding points to each cell, such that there are at
    ! least alpha_osc points per oscillatory wavelength and alpha_exp
    ! points per exponential wavelength.
    !
    ! Wavelengths are calculated based on a local dispersion analysis
    ! of the adibatic/Cowling wave equation, for inertial frequencies
    ! specified by omega

    allocate(rt, SOURCE=r_rot_t(ml, mp, op))

    ! At each point of x, determine the maximum absolute value of the
    ! real and imaginary parts of the local wavenumber

    n_x = SIZE(x)

    allocate(k_r_max(n_x))
    allocate(k_i_max(n_x))

    k_r_max(1) = 0._WP
    k_i_max(1) = 0._WP

    wavenumber_loop : do i = 2,n_x-1

       associate (V_g => ml%V_2(x(i))*x(i)**2/ml%Gamma_1(x(i)), As => ml%As(x(i)), &
                  U => ml%U(x(i)), c_1 => ml%c_1(x(i)))

         k_r_max(i) = 0._WP
         k_i_max(i) = 0._WP

         omega_loop : do j = 1, SIZE(omega)

            omega_c = rt%omega_c(x(i), omega(j))

            lambda = rt%lambda(x(i), omega(j))
            l_0 = rt%l_0(omega(j))
            
            ! Calculate the propagation discriminant gamma

            g_4 = -4._WP*V_g*c_1
            g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
            g_0 = -4._WP*lambda*As/c_1

            gamma = (g_4*omega_c**4 + g_2*omega_c**2 + g_0)/omega_c**2

            ! Update the wavenumber maxima

            if (gamma < 0._WP) then
               
               ! Propagation zone

               k_r_max(i) = MAX(k_r_max(i), ABS(0.5_WP*SQRT(-gamma))/x(i))
               k_i_max(i) = MAX(k_i_max(i), ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_0))/x(i))

            else

               ! Evanescent zone

               k_i_max(i) = MAX(k_i_max(i), ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_0 - SQRT(gamma)))/x(i), &
                                            ABS(0.5_WP*(As + V_g - U + 2._WP - 2._WP*l_0 + SQRT(gamma)))/x(i))

            end if

         end do omega_loop

       end associate

    end do wavenumber_loop

    k_r_max(n_x) = 0._WP
    k_i_max(n_x) = 0._WP

    ! Determine how many points to add to each cell

    allocate(dn(n_x-1))

    sample_loop : do i = 1,n_x-1

       ! Calculate the oscillatory and exponential phase change across
       ! the cell

       dphi_osc = MAX(k_r_max(i), k_r_max(i+1))*(x(i+1) - x(i))
       dphi_exp = MAX(k_i_max(i), k_i_max(i+1))*(x(i+1) - x(i))

       ! Set up dn

       dn(i) = MAX(FLOOR((alpha_osc*dphi_osc)/TWOPI), FLOOR((alpha_exp*dphi_exp)/TWOPI))

    end do sample_loop

    ! Perform the resampling

    call resample_(x, dn)

    ! Finish

    return

  end subroutine resample_dispersion_

!****

  subroutine resample_thermal_ (ml, mp, op, omega, alpha_thm, x)

    class(model_t), pointer, intent(in)  :: ml
    type(mode_par_t), intent(in)         :: mp
    type(osc_par_t), intent(in)          :: op
    real(WP), intent(in)                 :: omega(:)
    real(WP), intent(in)                 :: alpha_thm
    real(WP), allocatable, intent(inout) :: x(:)

    class(r_rot_t), allocatable :: rt
    integer                     :: n_x
    real(WP), allocatable       :: k_t_max(:)
    integer                     :: i
    integer                     :: j
    real(WP)                    :: omega_c
    integer, allocatable        :: dn(:)
    real(WP)                    :: dphi_thm

    $ASSERT(ALLOCATED(x),No input grid)

    ! Resample x by adding points to each cell, such that there are at
    ! least alpha_thm points per thermal length c_thm

    allocate(rt, SOURCE=r_rot_t(ml, mp, op))

    ! At each point of x, determine the maximum absolute value of the
    ! local thermal wavenumber

    n_x = SIZE(x)

    allocate(k_t_max(n_x))

    k_t_max(1) = 0._WP

    wavenumber_loop : do i = 2,n_x-1

       associate(V => ml%V_2(x(i))*x(i)**2, nabla => ml%nabla(x(i)), &
                c_rad => ml%c_rad(x(i)), c_thm => ml%c_thm(x(i)))
         
         k_t_max(i) = 0._WP

         omega_loop : do j = 1, SIZE(omega)

            omega_c = rt%omega_c(x(i), omega(j))

            k_t_max(i) = MAX(k_t_max(i), SQRT(ABS(V*nabla*omega_c*c_thm/c_rad))/x(i))

         end do omega_loop

       end associate

    end do wavenumber_loop

    k_t_max(n_x) = 0._WP

    ! Determine how many points to add to each cell

    allocate(dn(n_x-1))

    sample_loop : do i = 1,n_x-1

       ! Calculate the thermal phase change across the cell

       dphi_thm = MAX(k_t_max(i), k_t_max(i+1))*(x(i+1) - x(i))

       ! Set up dn

       dn(i) = FLOOR((alpha_thm*dphi_thm)/TWOPI)

    end do sample_loop

    ! Perform the resampling

    call resample_(x, dn)

    ! Finish

    return

  end subroutine resample_thermal_

!****

  subroutine resample_center_ (ml, mp, op, omega, n, x)

    class(model_t), pointer,  intent(in) :: ml
    type(mode_par_t), intent(in)         :: mp
    type(osc_par_t), intent(in)          :: op
    real(WP), intent(in)                 :: omega(:)
    integer, intent(in)                  :: n
    real(WP), allocatable, intent(inout) :: x(:)

    real(WP)             :: x_turn
    integer              :: j
    real(WP)             :: x_turn_j
    integer              :: i_turn
    integer              :: n_add
    integer              :: n_x
    integer, allocatable :: dn(:)

    $ASSERT(ALLOCATED(x),No input grid)

    ! Resample x by adding points to central cells, such that there
    ! are n_floor points covering the evanescent region at the
    ! center. Evanescence is determined based on a local dispersion
    ! analysis of the adibatic/Cowling wave equation, for frequencies
    ! in the interval [omega_min,omega_max]

    ! First, locate the innermost turning point

    x_turn = HUGE(0._WP)

    omega_loop : do j = 1, SIZE(omega)
       call find_x_turn(x, ml, mp, op, omega(j), x_turn_j)
       x_turn = MIN(x_turn, x_turn_j)
    end do omega_loop

    call locate(x, x_turn, i_turn)

    ! Determine how many points need to be added

    n_add = MAX(n-i_turn, 0)

    ! Determine how many points to add to each cell

    n_x = SIZE(x)

    allocate(dn(n_x-1))

    dn = 0

    if(i_turn >= 1 .AND. i_turn < n_x) then
       dn(:i_turn) = CEILING(n_add*(x(i_turn+1)/x_turn)/i_turn)
    endif
    
    ! Perform the resampling

    call resample_(x, dn)

    ! Finish

    return

  end subroutine resample_center_

!****

  subroutine resample_struct_ (ml, alpha_str, x)

    class(model_t), pointer, intent(in)  :: ml
    real(WP), intent(in)                 :: alpha_str
    real(WP), allocatable, intent(inout) :: x(:)

    integer :: dn(SIZE(x)-1)
    integer :: i

    $ASSERT(ALLOCATED(x),No input grid)

    ! Resample x by adding points to each cell, such that there are at
    ! least alpha_str points per dex change in the structure variables (V,
    ! As, Gamma_1, c_1, & U)

    ! Calculate the number of points to add to each cell

    dn = 0

    cell_loop : do i = 1, SIZE(x)-1

       dn(i) = dn(i) + FLOOR(alpha_str*dlog_(ml%V_2(x(i)), ml%V_2(x(i+1)))) + &
                       FLOOR(alpha_str*dlog_(ml%As(x(i)), ml%As(x(i+1)))) + &
                       FLOOR(alpha_str*dlog_(ml%Gamma_1(x(i)), ml%Gamma_1(x(i+1)))) + &
                       FLOOR(alpha_str*dlog_(ml%c_1(x(i)), ml%c_1(x(i+1)))) + &
                       FLOOR(alpha_str*dlog_(ml%U(x(i)), ml%U(x(i+1))))

    end do cell_loop

    ! Perform the resampling

    call resample_(x, dn)

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

  subroutine resample_uniform_ (n, x)

    integer, intent(in)                  :: n
    real(WP), allocatable, intent(inout) :: x(:)

    integer              :: n_x
    integer, allocatable :: dn(:)

    $ASSERT(ALLOCATED(x),No input grid)

    ! Resample x by adding n points to each cell

    n_x = SIZE(x)

    allocate(dn(n_x-1))

    dn = n

    ! Perform the resampling

    call resample_(x, dn)

    ! Finish

    return

  end subroutine resample_uniform_

!****

  subroutine find_x_turn (x, ml, mp, op, omega, x_turn)

    real(WP), intent(in)                :: x(:)
    class(model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)        :: mp
    type(osc_par_t), intent(in)         :: op
    real(WP), intent(in)                :: omega
    real(WP)                            :: x_turn

    type(gamma_func_t) :: gf
    integer            :: i

    ! Find the inner turning point at frequency omega

    x_turn = HUGE(0._WP)

    gf%ml => ml
    allocate(gf%rt, SOURCE=r_rot_t(ml, mp, op))

    gf%omega = omega

    turn_loop : do i = 1,SIZE(x)-1
       if(.NOT. (ml%is_zero(x(i)) .OR. ml%is_zero(x(i+1)))) then
          if(gf%eval(x(i)) > 0._WP .AND. gf%eval(x(i+1)) <= 0._WP) then
             x_turn = gf%root(x(i), x(i+1), 0._WP)
             exit turn_loop
          end if
       endif
    end do turn_loop

    ! Finish

    return

  end subroutine find_x_turn

!****

  function eval_gamma_func_ (this, z) result (gamma)

    class(gamma_func_t), intent(inout) :: this
    complex(WP), intent(in)            :: z
    complex(WP)                        :: gamma

    real(WP) :: x
    real(WP) :: g_4
    real(WP) :: g_2
    real(WP) :: g_0

    ! Calculate the propagation discriminant

    x = REAL(z)

    associate(V_g => this%ml%V_2(x)*x**2/this%ml%Gamma_1(x), As => this%ml%As(x), &
              U => this%ml%U(x), c_1 => this%ml%c_1(x), &
              lambda => this%rt%lambda(x, this%omega))

      g_4 = -4._WP*V_g*c_1
      g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
      g_0 = -4._WP*lambda*As/c_1

      gamma = (g_4*this%omega**4 + g_2*this%omega**2 + g_0)/this%omega**2

    end associate

    ! Finish

    return

  end function eval_gamma_func_

end module gyre_grid
