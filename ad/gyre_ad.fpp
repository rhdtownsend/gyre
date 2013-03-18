! Program  : gyre_ad
! Purpose  : adiabatic oscillation code
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

program gyre_ad

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  use core_hgroup
  use core_order

  use gyre_mech_coeffs
  use gyre_oscpar
  use gyre_ad_bvp
  use gyre_ad_search
  use gyre_mode
  use gyre_util
  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  integer                           :: unit
  real(WP), allocatable             :: x_mc(:)
  class(mech_coeffs_t), allocatable :: mc
  type(oscpar_t)                    :: op
  character(LEN=256)                :: ivp_solver_type
  type(ad_bvp_t)                    :: bp 
  real(WP), allocatable             :: omega(:)
  integer                           :: n_iter_max
  type(mode_t), allocatable         :: md(:)

  ! Initialize

  call init_parallel()

  if(MPI_RANK == 0) then

     write(OUTPUT_UNIT, '(A)') 'gyre_ad [hg]'
     write(OUTPUT_UNIT, '(A,2X,I0)') 'OpenMP Threads :', OMP_SIZE_MAX
     write(OUTPUT_UNIT, '(A,2X,I0)') 'MPI Processors :', MPI_SIZE
     
     call open_input(unit)

  endif

  call write_header('Initialization', '=')

  call init_coeffs(unit, x_mc, mc)
  call init_oscpar(unit, op)
  call init_numpar(unit, n_iter_max, ivp_solver_type)
  call init_scan(unit, mc, omega)
  call init_bvp(unit, x_mc, mc, op, omega, ivp_solver_type, bp)

  ! Find modes

  call ad_scan_search(bp, omega, n_iter_max, md)

  ! Write output
 
  if(MPI_RANK == 0) then
     call write_eigdata(unit, mc, op, md)
  endif

  ! Finish

  call final_parallel()

contains

  subroutine init_coeffs (unit, x_mc, mc)

    use gyre_mech_coeffs_evol
    use gyre_mech_coeffs_poly
    use gyre_mech_coeffs_hom
    use gyre_mesa_file
    use gyre_fgong_file
    use gyre_osc_file
    use gyre_b3_file

    integer, intent(in)                              :: unit
    real(WP), allocatable, intent(out)               :: x_mc(:)
    $if($GFORTRAN_PR56218)
    class(mech_coeffs_t), allocatable, intent(inout) :: mc
    $else
    class(mech_coeffs_t), allocatable, intent(out)   :: mc
    $endif

    character(LEN=256)          :: coeffs_type
    character(LEN=FILENAME_LEN) :: file
    real(WP)                    :: G
    real(WP)                    :: Gamma_1

    namelist /coeffs/ coeffs_type, file, G, Gamma_1

    ! Read structure coefficients parameters

    if(MPI_RANK == 0) then

       coeffs_type = ''

       file = ''

       G = G_GRAVITY
       Gamma_1 = 5._WP/3._WP

       rewind(unit)
       read(unit, NML=coeffs)

    endif

    $if($MPI)
    call bcast(coeffs_type, 0)
    $endif

    ! Allocate the mech_coeffs

    select case(coeffs_type)
    case('MESA','FGONG','B3')
       allocate(mech_coeffs_evol_t::mc)
    case('POLY')
       allocate(mech_coeffs_poly_t::mc)
    case('HOM')
       allocate(mech_coeffs_hom_t::mc)
    case default
       $ABORT(Invalid coeffs_type)
    end select

    ! Read the mech_coeffs

    if(MPI_RANK == 0) then

       select type (mc)
       type is (mech_coeffs_evol_t)
          select case(coeffs_type)
          case('MESA')
             call read_mesa_file(file, G, mc, x=x_mc)
          case('FGONG')
             call read_fgong_file(file, G, mc, x=x_mc)
          case('OSC')
             call read_osc_file(file, G, mc, x=x_mc)
          case('B3')
             call read_b3_file(file, G, mc, x=x_mc)
          end select
       type is (mech_coeffs_poly_t)
          call mc%read(file, x_mc)
       type is (mech_coeffs_hom_t)
          call mc%init(Gamma_1)
       end select

    endif

    $if($MPI)
    call mc%bcast(0)
    call bcast(x_mc, 0, alloc=.TRUE.)
    $endif

    ! Finish

    return

  end subroutine init_coeffs

!****

  subroutine init_oscpar (unit, op)

    integer, intent(in)         :: unit
    type(oscpar_t), intent(out) :: op

    integer            :: l
    character(LEN=256) :: outer_bound_type

    namelist /oscpar/ l, outer_bound_type

    ! Read oscillation parameters

    if(MPI_RANK == 0) then

       l = 0

       outer_bound_type = 'ZERO'

       rewind(unit)
       read(unit, NML=oscpar)

    endif

    $if($MPI)
    call bcast(l, 0)
    call bcast(outer_bound_type, 0)
    $endif

    ! Initialize the oscilaltion parameters

    call op%init(l, outer_bound_type)

    ! Finish

    return

  end subroutine init_oscpar

!****

  subroutine init_numpar (unit, n_iter_max, ivp_solver_type)

    integer, intent(in)           :: unit
    integer, intent(out)          :: n_iter_max
    character(LEN=*), intent(out) :: ivp_solver_type

    namelist /numpar/ n_iter_max, ivp_solver_type

    ! Read numerical parameters

    if(MPI_RANK == 0) then

       n_iter_max = 50

       ivp_solver_type = 'MAGNUS_GL2'

       rewind(unit)
       read(unit, NML=numpar)

    endif

    $if($MPI)
    call bcast(n_iter_max, 0)
    call bcast(ivp_solver_type, 0)
    $endif

    ! Finish

    return

  end subroutine init_numpar

!****

  subroutine init_scan (unit, mc, omega)

    integer, intent(in)                :: unit
    class(mech_coeffs_t), intent(in)   :: mc
    real(WP), allocatable, intent(out) :: omega(:)

    character(LEN=256) :: grid_type
    real(WP)           :: freq_min
    real(WP)           :: freq_max
    integer            :: n_freq
    character(LEN=256) :: freq_units
    real(WP)           :: omega_min
    real(WP)           :: omega_max
    integer            :: i

    namelist /scan/ grid_type, freq_min, freq_mid, freq_max, n_freq, freq_units, &
                    n_iter_max

    ! Read scan parameters

    if(MPI_RANK == 0) then

       rewind(unit)

       allocate(omega(0))

       read_loop : do 

          grid_type = 'LINEAR'

          freq_min = 1._WP
          freq_max = 10._WP
          n_freq = 10
          
          freq_units = 'NONE'

          read(unit, NML=scan, end=100)
          
          ! Set up the frequency grid

          omega_min = REAL(mc%conv_freq(CMPLX(freq_min, KIND=WP), freq_units, 'NONE'))
          omega_max = REAL(mc%conv_freq(CMPLX(freq_max, KIND=WP), freq_units, 'NONE'))
       
          select case(grid_type)
          case('LINEAR')
             omega = [omega,(((n_freq-i)*omega_min + (i-1)*omega_max)/(n_freq-1), i=1,n_freq)]
          case('INVERSE')
             omega = [omega,((n_freq-1)/((n_freq-i)/omega_min + (i-1)/omega_max), i=1,n_freq)]
          case default
             $ABORT(Invalid freq_grid)
          end select

       end do read_loop

100    continue

       ! Sort the frequencies

       omega = omega(sort_indices(omega))

    endif

    $if($MPI)
    call bcast(omega, 0, alloc=.TRUE.)
    $endif

    ! Finish

    return

  end subroutine init_scan

!****

  subroutine init_bvp (unit, x_mc, mc, op, omega, ivp_solver_type, bp)

    use gyre_ad_bound
    use gyre_ad_shooter
    use gyre_ad_jacobian

    integer, intent(in)                      :: unit
    real(WP), intent(in), allocatable        :: x_mc(:)
    class(mech_coeffs_t), intent(in), target :: mc
    type(oscpar_t), intent(in)               :: op
    real(WP), intent(in)                     :: omega(:)
    character(LEN=*), intent(in)             :: ivp_solver_type
    type(ad_bvp_t), intent(out)              :: bp

    type(ad_jacobian_t)   :: jc
    type(ad_bound_t)      :: bd
    character(LEN=256)    :: grid_type
    real(WP)              :: alpha_osc
    real(WP)              :: alpha_exp
    integer               :: n_center
    integer               :: n_floor
    real(WP)              :: s
    integer               :: n_grid
    integer               :: dn(SIZE(x_mc)-1)
    integer               :: i
    real(WP), allocatable :: x_sh(:)
    type(ad_shooter_t)    :: sh

    namelist /shoot_grid/ grid_type, alpha_osc, alpha_exp, &
         n_center, n_floor, s, n_grid

    namelist /recon_grid/ alpha_osc, alpha_exp, n_center, n_floor

    ! Initialize the Jacobian and boundary conditions

    call jc%init(mc, op)
    call bd%init(mc, op)

    ! Read shooting grid parameters

    if(MPI_RANK == 0) then

       grid_type = 'DISPERSION'

       alpha_osc = 0._WP
       alpha_exp = 0._WP
       
       n_center = 0
       n_floor = 0

       s = 100._WP
       n_grid = 100

       rewind(unit)
       read(unit, NML=shoot_grid)

    endif

    $if($MPI)
    call bcast(grid_type, 0)
    call bcast(alpha_osc, 0)
    call bcast(alpha_exp, 0)
    call bcast(n_center, 0)
    call bcast(n_floor, 0)
    call bcast(s, 0)
    call bcast(n_grid, 0)
    $endif

    ! Initialize the shooting grid

    select case(grid_type)
    case('GEOM')
       call build_geom_grid(s, n_grid, x_sh)
    case('LOG')
       call build_log_grid(s, n_grid, x_sh)
    case('INHERIT')
       $ASSERT(ALLOCATED(x_mc),No input grid)
       x_sh = x_mc
    case('DISPERSION')
       $ASSERT(ALLOCATED(x_mc),No input grid)
       dn = 0
       do i = 1,SIZE(omega)
          call plan_dispersion_grid(x_mc, mc, CMPLX(omega(i), KIND=WP), op, alpha_osc, alpha_exp, n_center, n_floor, dn)
       enddo
       call build_oversamp_grid(x_mc, dn, x_sh)
    case default
       $ABORT(Invalid grid_type)
    end select

    ! Read recon grid parameters

    if(MPI_RANK == 0) then

       alpha_osc = 0._WP
       alpha_exp = 0._WP
       
       n_center = 0
       n_floor = 0

       rewind(unit)
       read(unit, NML=recon_grid)

    endif

    $if($MPI)
    call bcast(alpha_osc, 0)
    call bcast(alpha_exp, 0)
    call bcast(n_center, 0)
    call bcast(n_floor, 0)
    $endif

    ! Initialize the shooter

    call sh%init(mc, op, jc, x_sh, alpha_osc, alpha_exp, n_center, n_floor, ivp_solver_type)

    ! Initialize the bvp

    call bp%init(sh, bd)

    ! Finish

    return

  end subroutine init_bvp

!****

  subroutine write_eigdata (unit, mc, op, md)

    integer, intent(in)              :: unit
    class(mech_coeffs_t), intent(in) :: mc
    class(oscpar_t), intent(in)      :: op
    type(mode_t), intent(in)         :: md(:)

    character(LEN=256)               :: freq_units
    character(LEN=FILENAME_LEN)      :: eigval_file
    character(LEN=FILENAME_LEN)      :: eigfunc_prefix
    integer                          :: i
    real(WP)                         :: E(SIZE(md))
    complex(WP)                      :: freq(SIZE(md))
    type(hgroup_t)                   :: hg
    character(LEN=FILENAME_LEN)      :: eigfunc_file

    namelist /output/ freq_units, eigval_file, eigfunc_prefix

    ! Read output parameters

    freq_units = 'NONE'

    eigval_file = ''
    eigfunc_prefix = ''

    rewind(unit)
    read(unit, NML=output)

    ! Calculate inertias

    do i = 1,SIZE(md)
       E(i) = inertia(mc, op, md(i))
    end do

    ! Write eigenvalues

    freq_loop : do i = 1,SIZE(md)
       freq(i) = mc%conv_freq(REAL(md(i)%omega), 'NONE', freq_units)
    end do freq_loop

    if(eigval_file /= '') then

       call hg%init(eigval_file, CREATE_FILE)
          
       call write_attr(hg, 'n_md', SIZE(md))

       call write_attr(hg, 'n', bp%n)
       call write_attr(hg, 'n_e', bp%n_e)

       call write_attr(hg, 'l', op%l)

       call write_dset(hg, 'n_p', md%n_p)
       call write_dset(hg, 'n_g', md%n_g)

       call write_dset(hg, 'freq', freq)
       call write_attr(hg, 'freq_units', freq_units)

       call write_dset(hg, 'E', E)

       call hg%final()

    end if

    ! Write eigenfunctions

    if(eigfunc_prefix /= '') then

       mode_loop : do i = 1,SIZE(md)

          write(eigfunc_file, 100) TRIM(eigfunc_prefix), i, '.h5'
100       format(A,I4.4,A)

          call hg%init(eigfunc_file, CREATE_FILE)

          call write_attr(hg, 'n', bp%n)
          call write_attr(hg, 'n_e', bp%n_e)

          call write_attr(hg, 'l', op%l)
          call write_attr(hg, 'lambda_0', op%lambda_0)

          call write_attr(hg, 'n_p', md(i)%n_p)
          call write_attr(hg, 'n_g', md(i)%n_g)

          call write_attr(hg, 'freq', freq(i))
          call write_attr(hg, 'freq_units', freq_units)

          call write_dset(hg, 'x', md(i)%x)
          call write_dset(hg, 'y', md(i)%y)

          call write_dset(hg, 'dE_dx', kinetic(mc, op, md(i)))

          call hg%final()

       end do mode_loop

    end if

    ! Finish

    return

  end subroutine write_eigdata

!*****

  function inertia (mc, op, md) result (E)

    class(mech_coeffs_t), intent(in) :: mc
    class(oscpar_t), intent(in)      :: op
    class(mode_t), intent(in)        :: md
    real(WP)                         :: E

    real(WP)    :: dE_dx(md%n)
    complex(WP) :: xi_r
    complex(WP) :: xi_h
    real(WP)    :: E_norm

    ! Calculate the normalized inertia

    dE_dx = kinetic(mc, op, md)

    E = SUM(0.5_WP*(dE_dx(2:) + dE_dx(:md%n-1))*(md%x(2:) - md%x(:md%n-1)))

    if(op%l == 0) then
       xi_r = md%y(1,md%n)
       xi_h = 0._WP
    else
       xi_r = md%y(1,md%n)
       xi_h = md%y(2,md%n)/md%omega**2
    endif

    E_norm = ABS(xi_r)**2 + op%l*(op%l+1)*ABS(xi_h)**2
    $ASSERT(E_norm /= 0._WP,E_norm is zero)    
    
    E = E/E_norm

    ! Finish

    return

  end function inertia

!*****

  function kinetic (mc, op, md) result (dE_dx)

    class(mech_coeffs_t), intent(in) :: mc
    class(oscpar_t), intent(in)      :: op
    class(mode_t), intent(in)        :: md
    real(WP)                         :: dE_dx(md%n)

    integer     :: i
    real(WP)    :: U
    real(WP)    :: c_1
    complex(WP) :: xi_r
    complex(WP) :: xi_h

    ! Calculate the kinetic energy density

    !$OMP PARALLEL DO PRIVATE (U, c_1, xi_r, xi_h)
    do i = 1,md%n

       U = mc%U(md%x(i))
       c_1 = mc%c_1(md%x(i))

       if(op%l == 0) then
          xi_r = md%y(1,i)
          xi_h = 0._WP
       else
          xi_r = md%y(1,i)
          xi_h = md%y(2,i)/(c_1*md%omega**2)
       endif

       if(md%x(i) > 0._WP) then
          xi_r = xi_r*md%x(i)**(op%lambda_0+1._WP)
          xi_h = xi_h*md%x(i)**(op%lambda_0+1._WP)
       else
          if(op%lambda_0 /= -1._WP) then
             xi_r = 0._WP
             xi_h = 0._WP
          endif
       endif

       dE_dx(i) = (ABS(xi_r)**2 + op%l*(op%l+1)*ABS(xi_h)**2)*U*md%x(i)**2/c_1

    end do

    ! Finish

    return

  end function kinetic

end program gyre_ad
