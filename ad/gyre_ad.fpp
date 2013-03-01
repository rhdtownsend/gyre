! Program  : gyre_ad
! Purpose  : Adiabatic oscillations

$include 'core.inc'

program gyre_ad

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  use core_hgroup

  use gyre_input
  use gyre_mech_coeffs
  use gyre_mech_coeffs_evol
  use gyre_mech_coeffs_poly
  use gyre_mech_coeffs_hom
  use gyre_ad_bvp
  use gyre_ad_api
  use gyre_ad_bound
  use gyre_ad_shooter
  use gyre_ad_jacobian
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
  type(ad_bvp_t)                    :: bp
  real(WP), allocatable             :: omega(:)
  integer                           :: n_iter_max
  type(mode_t), allocatable         :: md(:)
  integer                           :: i
  real(WP), allocatable             :: E(:)

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
  call init_bvp(unit, x_mc, mc, bp)
  call init_scan(unit, mc, omega, n_iter_max)

  ! Find modes

  call find_ad_modes(bp, omega, n_iter_max, md)

  ! Calculate inertias

  allocate(E(SIZE(md)))

  do i = 1,SIZE(md)
!     E(i) = inertia(md(i), bp)
     E(i) = 0._WP
  end do

  ! Write output
 
  if(MPI_RANK == 0) then
     call write_eigdata(unit, mc, md)
  endif

  ! Finish

  call final_parallel()

contains

  subroutine init_coeffs (unit, x_mc, mc)

    integer, intent(in)                              :: unit
    real(WP), allocatable, intent(out)               :: x_mc(:)
    $if($GFORTRAN_PR_56218)
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
          call mc%read(coeffs_type, file, G, x_mc)
       type is (mech_coeffs_poly_t)
          call mc%read(file, x_mc)
       type is (mech_coeffs_hom_t)
          call mc%init(Gamma_1)
       end select

    endif

    $if($MPI)
    call mc%bcast(0)
    $endif

    ! Finish

    return

  end subroutine init_coeffs

!****

  subroutine init_bvp (unit, x_mc, mc, bp)

    integer, intent(in)                      :: unit
    real(WP), intent(in), allocatable        :: x_mc(:)
    class(mech_coeffs_t), intent(in), target :: mc
    type(ad_bvp_t), intent(out)              :: bp

    integer               :: l
    character(LEN=256)    :: outer_bound
    character(LEN=256)    :: solver
    real(WP)              :: lambda_0
    type(ad_jacobian_t)   :: jc
    type(ad_bound_t)      :: bd
    character(LEN=256)    :: grid_type
    real(WP)              :: freq
    character(LEN=256)    :: freq_units
    real(WP)              :: alpha_osc
    real(WP)              :: alpha_exp
    integer               :: n_center
    integer               :: n_floor
    real(WP)              :: s
    integer               :: n_grid
    integer               :: dn(SIZE(x_mc)-1)
    complex(WP)           :: omega
    real(WP), allocatable :: x_sh(:)
    type(ad_shooter_t)    :: sh

    namelist /eqns/ l, outer_bound, solver

    namelist /shoot_grid/ grid_type, freq, freq_units, alpha_osc, alpha_exp, &
         n_center, n_floor, s, n_grid

    namelist /recon_grid/ alpha_osc, alpha_exp, n_center, n_floor

    ! Read equation parameters

    if(MPI_RANK == 0) then

       l = 0

       outer_bound = 'ZERO'
       solver = 'MAGNUS_GL2'

       rewind(unit)
       read(unit, NML=eqns)

    endif

    $if($MPI)
    call bcast(l, 0)
    call bcast(outer_bound, 0)
    call bcast(solver, 0)
    $endif

    ! Initialize the Jacobian and boundary conditions

    if(l == 0) then
       lambda_0 = 0._WP
    else
       lambda_0 = l - 2._WP
    endif

    call jc%init(mc, lambda_0, l)
    call bd%init(mc, lambda_0, l, outer_bound)

    ! Read shooting grid parameters

    if(MPI_RANK == 0) then

       grid_type = 'DISPERSION'

       freq = 1._WP
       freq_units = 'NONE'

       alpha_osc = 0._WP
       alpha_exp = 0._WP
       
       n_center = 0
       n_floor = 0

       s = 0._WP
       n_grid = 0

       rewind(unit)
       read(unit, NML=shoot_grid)

    endif

    $if($MPI)
    call bcast(grid_type, 0)
    call bcast(freq, 0)
    call bcast(freq_units, 0)
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
    case('DISPERSION')
       $ASSERT(ALLOCATED(x_mc),No input grid)
       omega = mc%conv_freq(freq, freq_units, 'NONE')
       call plan_dispersion_grid(x_mc, mc, omega, l, alpha_osc, alpha_exp, n_center, n_floor, dn)
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

    call sh%init(mc, jc, x_sh, l, alpha_osc, alpha_exp, n_center, n_floor, solver)

    ! Initialize the bvp

    call bp%init(sh, bd)

    ! Finish

    return

  end subroutine init_bvp

!****

  subroutine init_scan (unit, mc, omega, n_iter_max)

    integer, intent(in)                :: unit
    class(mech_coeffs_t), intent(in)   :: mc
    real(WP), allocatable, intent(out) :: omega(:)
    integer, intent(out)               :: n_iter_max

    character(LEN=256) :: grid_type
    real(WP)           :: freq_min
    real(WP)           :: freq_mid
    real(WP)           :: freq_max
    integer            :: n_freq
    character(LEN=256) :: freq_units
    real(WP)           :: omega_min
    real(WP)           :: omega_mid
    real(WP)           :: omega_max
    integer            :: i

    namelist /scan/ grid_type, freq_min, freq_mid, freq_max, n_freq, freq_units, &
                    n_iter_max

    ! Read scan parameters

    if(MPI_RANK == 0) then

       grid_type = 'LINEAR'

       freq_min = 1._WP
       freq_mid = SQRT(10._WP)
       freq_max = 10._WP
       n_freq = 10
    
       freq_units = 'NONE'

       n_iter_max = 50

       rewind(unit)
       read(unit, NML=scan)

    endif

    ! Set up the frequency grid

    if(MPI_RANK == 0) then

       omega_min = mc%conv_freq(freq_min, freq_units, 'NONE')
       omega_mid = mc%conv_freq(freq_mid, freq_units, 'NONE')
       omega_max = mc%conv_freq(freq_max, freq_units, 'NONE')

       select case(grid_type)
       case('LINEAR')
          omega = [(((n_freq-i)*omega_min + (i-1)*omega_max)/(n_freq-1), i=1,n_freq)]
       case('INVERSE')
          omega = [((n_freq-1)/((n_freq-i)/omega_min + (i-1)/omega_max), i=1,n_freq)]
       case('MIXED')
          omega = [((n_freq-1)/((n_freq-i)/omega_min + (i-1)/omega_mid), i=1,n_freq), &
               (((n_freq-i)*omega_mid + (i-1)*omega_max)/(n_freq-1), i=1,n_freq)]
       case default
          $ABORT(Invalid freq_grid)
       end select

    endif

    $if($MPI)
    call bcast(omega, 0, alloc=.TRUE.)
    call bcast(n_iter_max, 0)
    $endif

    ! Finish

    return

  end subroutine init_scan

!****

  subroutine write_eigdata (unit, mc, md)

    integer, intent(in)              :: unit
    class(mech_coeffs_t), intent(in) :: mc
    type(mode_t), intent(in)         :: md(:)

    character(LEN=256)               :: freq_units
    character(LEN=FILENAME_LEN)      :: eigval_file
    character(LEN=FILENAME_LEN)      :: eigfunc_prefix
    integer                          :: i
    complex(WP)                      :: freq(SIZE(md))
    type(hgroup_t)                   :: hg
    character(LEN=FILENAME_LEN)      :: eigfunc_file

    namelist /output/ freq_units, eigval_file, eigfunc_prefix

    ! Set defaults

    freq_units = 'NONE'

    eigval_file = ''
    eigfunc_prefix = ''

    ! Read parameters

    rewind(unit)
    read(unit, NML=output)

    ! Write eigenvalues

    freq_loop : do i = 1,SIZE(md)
       freq(i) = mc%conv_freq(REAL(md(i)%omega), 'NONE', freq_units)
    end do freq_loop

    if(eigval_file /= '') then

       call hg%init(eigval_file, CREATE_FILE)
          
       call write_attr(hg, 'n_md', SIZE(md))

       call write_attr(hg, 'n', bp%n)
       call write_attr(hg, 'n_e', bp%n_e)

!       call write_attr(hg, 'l', bp%l)

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

!          call write_attr(hg, 'l', bp%l)
!          call write_attr(hg, 'lambda_0', bp%lambda_0)

          call write_attr(hg, 'n_p', md(i)%n_p)
          call write_attr(hg, 'n_g', md(i)%n_g)

          call write_attr(hg, 'freq', freq(i))
          call write_attr(hg, 'freq_units', freq_units)

          call write_dset(hg, 'x', md(i)%x)
          call write_dset(hg, 'y', md(i)%y)

          call hg%final()

       end do mode_loop

    end if

    ! Finish

    return

  end subroutine write_eigdata

end program gyre_ad
