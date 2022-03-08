! Program  : scan_overlap_res
! Purpose  : scan through resonances using overlap method
!
! Copyright 2021 Rich Townsend & The GYRE Team
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

program scan_overlap_res
   
   ! Uses
   
   use core_kinds
   use core_system
   use core_parallel
   use core_order
   use core_memory
   use core_hgroup

   use gyre_constants
   use gyre_func
   use gyre_interp
   use gyre_math
   use gyre_mode_par
   use gyre_tide_util
   use gyre_util

   use ISO_FORTRAN_ENV

   ! No implicit typing

   implicit none

   ! Derived-type definitions

   type orbit_par_t
      real(WP) :: Omega_orb = 0._WP
      real(WP) :: Omega_rot = 0._WP
      real(WP) :: q = 1._WP
      real(WP) :: e = 0._WP
      integer  :: k = 0
   end type orbit_par_t

   type overlap_par_t
      character(32)  :: scan_var = ''
      integer        :: scan_n_pg_min = -1
      integer        :: scan_n_pg_max = -1
      integer        :: sum_n_pg_min = -1
      integer        :: sum_n_pg_max = -1
      real(WP)       :: del_gam  = 0.1_WP
      real(WP)       :: del_eps = 1.01_WP
      real(WP)       :: gamma_min = 0._WP
      logical        :: exclude_center = .FALSE.
      character(256) :: summary_file = ''
      character(256) :: detail_template = ''
      character(256) :: output_file = ''
   end type overlap_par_t

   type summary_data_t
      complex(WP), allocatable :: omega(:)
      integer, allocatable     :: n_pg(:)
      integer, allocatable     :: j(:)
      integer                  :: n
   end type summary_data_t

   type detail_data_t
      complex(WP)              :: omega
      real(WP), allocatable    :: x(:)
      complex(WP), allocatable :: xi_r(:)
      complex(WP), allocatable :: xi_h(:)
      complex(WP), allocatable :: eul_phi(:)
      real(WP)                 :: E
      integer                  :: n_pg
      integer                  :: j
      integer                  :: n
   end type detail_data_t

   type response_t
      complex(WP)              :: omega
      real(WP), allocatable    :: sigma(:)
      real(WP), allocatable    :: Omega_orb(:)
      real(WP), allocatable    :: Omega_rot(:)
      complex(WP), allocatable :: xi_r(:)
      complex(WP), allocatable :: xi_h(:)
      complex(WP), allocatable :: eul_phi(:)
      complex(WP), allocatable :: F(:)
      real(WP), allocatable    :: tau(:)
      complex(WP), allocatable :: xi_r_res(:)
      complex(WP), allocatable :: xi_h_res(:)
      complex(WP), allocatable :: eul_phi_res(:)
      complex(WP), allocatable :: F_res(:)
      real(WP), allocatable    :: tau_res(:)
      integer                  :: n
      integer                  :: n_pg
   end type response_t

   ! Variables

   character(:), allocatable        :: filename
   integer                          :: unit
   type(mode_par_t), allocatable    :: md_p(:)
   type(orbit_par_t)                :: or_p
   type(overlap_par_t)              :: ov_p
   type(summary_data_t)             :: sd
   type(detail_data_t), allocatable :: dd(:)
   integer                          :: n_res
   type(response_t), allocatable    :: rs(:)
   integer                          :: i
   integer                          :: n_pg
   type(hgroup_t)                   :: hg
   type(hgroup_t)                   :: hg_rs

   ! Read command-line arguments

   $ASSERT(n_arg() == 1,Syntax: scan_overlap_res <filename>)

   call get_arg(1, filename)

   ! Initialize

   call init_parallel()
   call init_math()

   ! Read the namelist file

   open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

   call read_constants(unit)

   call read_mode_par(unit, md_p)
   call read_orbit_par(unit, or_p)
   call read_overlap_par(unit, ov_p)

   $ASSERT(SIZE(md_p) == 1,Only a single &mode namelist allowed)

   ! Read data from the summary file

   call read_summary_data(ov_p, sd)

   ! Read data from the detail files

   call read_detail_data(ov_p, md_p(1), sd, dd)

   ! Allocate space for results

   n_res= ov_p%scan_n_pg_max - ov_p%scan_n_pg_min + 1

   allocate(rs(n_res))

   ! Now loop through the resonance points

   do i = 1, n_res

      ! Determine the resonance we are scanning

      n_pg = ov_p%scan_n_pg_min + i - 1

      ! Create the response

      call create_response(sd, dd, n_pg, md_p(1), or_p, ov_p, rs(i))

   end do

   ! Write out results

   hg = hgroup_t(ov_p%output_file, CREATE_FILE)

   call write_attr(hg, 'n_res', n_res)

   do i = 1, n_res

      hg_rs = hgroup_t(hg, rs_name(i))
      call write_response(hg_rs, rs(i))
      call hg_rs%final()

   end do

   call hg%final()

   ! Finish

   close(unit)

   call final_parallel()

contains

   subroutine read_orbit_par(unit, or_p)

      integer, intent(in)            :: unit
      type(orbit_par_t), intent(out) :: or_p

      integer  :: n_or_p
      real(WP) :: Omega_orb
      real(WP) :: Omega_rot
      real(WP) :: q
      real(WP) :: e
      integer  :: k

      namelist /orbit/ Omega_orb, Omega_rot, q, e, k

      ! Count the number of orbit namelists

      rewind(unit)

      n_or_p = 0

      count_loop : do
         read(unit, NML=orbit, END=100)
         n_or_p = n_or_p + 1
      end do count_loop

100   continue

      $ASSERT(n_or_p == 1,Must be exacly one &orbit namelist)

      ! Read orbit parameters

      rewind(unit)

      ! Set default values

      or_p = orbit_par_t()

      Omega_orb = or_p%Omega_orb
      Omega_rot = or_p%Omega_rot
      q = or_p%q
      e = or_p%e
      k = or_p%k

      ! Read the namelist

      read(unit, NML=orbit)

      ! Store read values

      or_p%Omega_orb = Omega_orb
      or_p%Omega_rot = Omega_rot
      or_p%q = q
      or_p%e = e
      or_p%k = k

      ! Finish

      return

   end subroutine read_orbit_par

   !****

   subroutine read_overlap_par(unit, ov_p)

      integer, intent(in)              :: unit
      type(overlap_par_t), intent(out) :: ov_p

      integer        :: n_ov_p
      character(32)  :: scan_var
      integer        :: scan_n_pg_min
      integer        :: scan_n_pg_max
      integer        :: sum_n_pg_min
      integer        :: sum_n_pg_max
      real(WP)       :: del_gam
      real(WP)       :: del_eps
      real(WP)       :: gamma_min
      logical        :: exclude_center
      character(256) :: summary_file
      character(256) :: detail_template
      character(256) :: output_file

      namelist /overlap/ scan_var, scan_n_pg_min, scan_n_pg_max, sum_n_pg_min, sum_n_pg_max, &
           del_gam, del_eps, gamma_min, exclude_center, summary_file, detail_template, output_file

      ! Count the number of overlap namelists

      rewind(unit)

      n_ov_p = 0

      count_loop : do
         read(unit, NML=overlap, END=100)
         n_ov_p = n_ov_p + 1
      end do count_loop

100   continue

      $ASSERT(n_ov_p == 1,Must be exacly one &overlap namelist)

      ! Read overlap parameters

      rewind(unit)

      ! Set default values

      ov_p = overlap_par_t()

      scan_var = ov_p%scan_var
      scan_n_pg_min = ov_p%scan_n_pg_min
      scan_n_pg_max = ov_p%scan_n_pg_max
      sum_n_pg_min = ov_p%sum_n_pg_min
      sum_n_pg_max = ov_p%sum_n_pg_max
      del_gam = ov_p%del_gam
      del_eps = ov_p%del_eps
      gamma_min = ov_p%gamma_min
      exclude_center = ov_p%exclude_center
      summary_file = ov_p%summary_file
      detail_template = ov_p%detail_template
      output_file = ov_p%output_file

      ! Read the namelist

      read(unit, NML=overlap)

      ! Store read values

      ov_p%scan_var = scan_var
      ov_p%scan_n_pg_min = scan_n_pg_min
      ov_p%scan_n_pg_max = scan_n_pg_max
      ov_p%sum_n_pg_min = sum_n_pg_min
      ov_p%sum_n_pg_max = sum_n_pg_max
      ov_p%del_gam = del_gam 
      ov_p%del_eps = del_eps
      ov_p%gamma_min = gamma_min
      ov_p%exclude_center = exclude_center
      ov_p%summary_file = summary_file
      ov_p%detail_template = detail_template
      ov_p%output_file = output_file

      ! Finish

      return

   end subroutine read_overlap_par

   !****

   subroutine read_summary_data(ov_p, sd)

      type(overlap_par_t), intent(in)   :: ov_p
      type(summary_data_t), intent(out) :: sd

      type(hgroup_t)           :: hg
      integer, allocatable     :: n_pg(:)
      integer, allocatable     :: j(:)
      complex(WP), allocatable :: omega(:)
      logical, allocatable     :: mask(:)
      integer, allocatable     :: i_srt(:)
      integer                  :: n
      integer                  :: i
      integer, allocatable     :: n_pg_dd(:)
      complex(WP), allocatable :: omega_dd(:)
      type(r_interp_t)         :: omega_in
      integer                  :: k

      ! Read the summary data

      print *,'Read summary:', ov_p%summary_file

      hg = hgroup_t(ov_p%summary_file, OPEN_FILE)

      call read_dset_alloc(hg, 'n_pg_ad', n_pg)
      call read_dset_alloc(hg, 'j', j)
      call read_dset_alloc(hg, 'omega', omega)

      call hg%final()

      ! Perform a bit of quality-control on the data

      ! Limit data to the range indicated by sum_n_pg_[min/max]

      mask = n_pg >= ov_p%sum_n_pg_min .AND. n_pg <= ov_p%sum_n_pg_max

      n_pg = PACK(n_pg, mask)
      j = PACK(j, mask)
      omega = PACK(omega, mask)

      n = SIZE(n_pg)

      print *,'  ...limited to', n

      ! Sort

      i_srt = sort_indices(n_pg)

      n_pg =  n_pg(i_srt)
      j = j(i_srt)
      omega = omega(i_srt)

      print *,'  ...sorted'

      ! Set up de-duplicated arrays

      deallocate(mask)
      allocate(mask(n))

      do i = 1, n
         mask(i) = COUNT(n_pg == n_pg(i)) == 1
      enddo

      n_pg_dd = PACK(n_pg, mask)
      omega_dd = PACK(omega, mask)

      print *,'  ...created dup-free', SIZE(n_pg_dd)

      ! Set up an interpolating function for the de-duped data

      omega_in = r_interp_t(REAL(n_pg_dd, KIND=WP), omega_dd%re, 'MONO')

      print *,'  ...created fit spline', MINVAL(n_pg_dd), MAXVAL(n_pg_dd)

      ! Now populate sd using data arrays with duplicates resolved
      ! using the interpolating function

      i_srt = unique_indices(n_pg)

      sd%n_pg = n_pg(i_srt)

      sd%n = SIZE(sd%n_pg)

      allocate(sd%j(sd%n))
      allocate(sd%omega(sd%n))

      do i = 1, sd%n

         if (COUNT(n_pg == sd%n_pg(i)) > 1) then
            k = MINLOC(ABS(omega_in%f(REAL(sd%n_pg(i), KIND=WP)) - omega%re), DIM=1)
         else
            k = FINDLOC(n_pg, sd%n_pg(i), DIM=1)
         endif

         sd%j(i) = j(k)
         sd%omega(i) = omega(k)

      end do

      print *,'  ...populated with dup res', sd%n

      ! Check for missing modes

      do i = 2, sd%n
         if (sd%n_pg(i) /= sd%n_pg(i-1)+1) then
            print *,'n_pg jumps from ', sd%n_pg(i-1), ' to ', sd%n_pg(i)
            $ABORT(Missing modes)
         endif
      enddo

      ! Finish

      return

   end subroutine read_summary_data

   !****

   subroutine read_detail_data(ov_p, md_p, sd, dd)

      type(overlap_par_t), intent(in)               :: ov_p
      type(mode_par_t), intent(in)                  :: md_p
      type(summary_data_t), intent(in)              :: sd
      type(detail_data_t), allocatable, intent(out) :: dd(:)

      integer                   :: n
      integer                   :: i
      integer                   :: n_pg
      character(:), allocatable :: detail_file
      type(hgroup_t)            :: hg

      ! Read the detail data

      allocate(dd(sd%n))

      do i = 1, sd%n

         ! Set up the filename

         detail_file = ov_p%detail_template

         detail_file = subst(detail_file, '%J', sd%j(i), '(I5.5)')
         detail_file = subst(detail_file, '%L', md_p%l, '(I3.3)')
         detail_file = subst(detail_file, '%M', md_p%m, '(SP,I3.2)')
         detail_file = subst(detail_file, '%j', sd%j(i), '(I0)')
         detail_file = subst(detail_file, '%l', md_p%l, '(I0)')
         detail_file = subst(detail_file, '%m', md_p%m, '(SP,I0)')

         ! Read the data

         print *,'Read detail:', detail_file

         hg = hgroup_t(detail_file, OPEN_FILE)

         call read_attr(hg, 'omega', dd(i)%omega)
         call read_attr(hg, 'n_pg', dd(i)%n_pg)
         call read_attr(hg, 'E', dd(i)%E)

         call read_dset_alloc(hg, 'x', dd(i)%x)
         call read_dset_alloc(hg, 'xi_r', dd(i)%xi_r)
         call read_dset_alloc(hg, 'xi_h', dd(i)%xi_h)
         call read_dset_alloc(hg, 'eul_phi', dd(i)%eul_phi)

         call hg%final()

         dd(i)%n = SIZE(dd(i)%x)

      end do

      $ASSERT(ALL(dd%n == dd(1)%n),Inconsistent dimensions)

      ! Finish

      return

   end subroutine read_detail_data

   !****

   subroutine create_response(sd, dd, n_pg, md_p, or_p, ov_p, rs)

      type(summary_data_t), intent(in) :: sd
      type(detail_data_t), intent(in)  :: dd(:)
      integer, intent(in)              :: n_pg
      type(mode_par_t), intent(in)     :: md_p
      type(orbit_par_t), intent(in)    :: or_p
      type(overlap_par_t), intent(in)  :: ov_p
      type(response_t), intent(out)    :: rs

      integer  :: i_sd
      integer  :: i_sd_m
      integer  :: i_sd_p
      real(WP) :: sigma_min
      real(WP) :: sigma_max
      integer  :: d
      integer  :: n
      real(WP) :: gamma
      real(WP) :: sigma
      real(WP) :: dsigma
      integer  :: j

      ! Find where in the summary data this and adjacent resonances
      ! are located

      i_sd = FINDLOC(sd%n_pg, n_pg, DIM=1)
      i_sd_m = FINDLOC(sd%n_pg, n_pg-1, DIM=1)
      i_sd_p = FINDLOC(sd%n_pg, n_pg+1, DIM=1)

      rs%omega = sd%omega(i_sd)
      rs%n_pg = n_pg

      ! Set up co-rotating frequency bounds

      sigma_min = 0.5_WP*(sd%omega(i_sd_m)%re + sd%omega(i_sd)%re)
      sigma_max = 0.5_WP*(sd%omega(i_sd)%re + sd%omega(i_sd_p)%re)

      ! Set up the co-rotating frequency grid

      d = 1024
      n = 0

      allocate(rs%sigma(d))

      if (ov_p%exclude_center) then
         n = 0
      else
         n = 1
         rs%sigma(1) = sd%omega(i_sd)%re
      endif

      gamma = MAX(ABS(sd%omega(i_sd)%im), ov_p%gamma_min)

      sigma = sd%omega(i_sd)%re
      dsigma = ov_p%del_gam*gamma

      do

         n = n + 1

         if (n > d) then
            d = 2*d
            call reallocate(rs%sigma, [d])
         endif

         sigma = sigma - dsigma

         if (sigma < sigma_min) then
            rs%sigma(n) = sigma_min
            exit
         else
            rs%sigma(n) = sigma
         endif

         dsigma = dsigma*ov_p%del_eps

      end do

      sigma = sd%omega(i_sd)%re
      dsigma = ov_p%del_gam*gamma

      do

         n = n + 1

         if (n > d) then
            d = 2*d
            call reallocate(rs%sigma, [d])
         endif

         sigma = sigma + dsigma

         if (sigma > sigma_max) then
            rs%sigma(n) = sigma_max
            exit
         else
            rs%sigma(n) = sigma
         endif

         dsigma = dsigma*ov_p%del_eps

      end do

      rs%sigma = rs%sigma(sort_indices(rs%sigma(:n)))

      rs%n = n

      print *,'scanning:', n_pg, rs%omega, rs%n

      ! Allocate arrays

      allocate(rs%Omega_orb(n))
      allocate(rs%Omega_rot(n))

      allocate(rs%xi_r(n))
      allocate(rs%xi_h(n))
      allocate(rs%eul_phi(n))
      allocate(rs%F(n))
      allocate(rs%tau(n))
 
      allocate(rs%xi_r_res(n))
      allocate(rs%xi_h_res(n))
      allocate(rs%eul_phi_res(n))
      allocate(rs%F_res(n))
      allocate(rs%tau_res(n))

      ! Set up Omega_orb and Omerga_rot

      select case (ov_p%scan_var)

      case ('OMEGA_ORB')

         $ASSERT(or_p%k /= 0,Cannot scan over Omega_orb when k == 0)

         rs%Omega_orb = (rs%sigma + md_p%m*or_p%Omega_rot)/or_p%k
         rs%Omega_rot = or_p%Omega_rot

      case ('OMEGA_ROT')

         $ASSERT(md_p%m /= 0,Cannot scan over Omega_rot when m == 0)

         rs%Omega_orb = or_p%Omega_orb
         rs%Omega_rot = (or_p%k*or_p%Omega_orb - rs%sigma)/md_p%m

      case default

         $ABORT(Invalid scan_var)

      end select

      ! Now evaluate the response data

      !$OMP PARALLEL DO
      do j = 1, n
         call eval_response(sd, dd, rs%Omega_orb(j), rs%Omega_rot(j), md_p, or_p, ov_p, j, rs)
      end do

      ! Finish

      return

   end subroutine create_response

   !****

   subroutine eval_response(sd, dd, Omega_orb, Omega_rot, md_p, or_p, ov_p, j, rs)

      type(summary_data_t), intent(in) :: sd
      type(detail_data_t), intent(in)  :: dd(:)
      real(WP), intent(in)             :: Omega_orb
      real(WP), intent(in)             :: Omega_rot
      type(mode_par_t), intent(in)     :: md_p
      type(orbit_par_t), intent(in)    :: or_p
      type(overlap_par_t), intent(in)  :: ov_p
      integer, intent(in)              :: j
      type(response_t), intent(inout)  :: rs

      integer     :: n
      integer     :: i
      real(WP)    :: omega
      real(WP)    :: gamma
      real(WP)    :: R_a
      real(WP)    :: eps
      real(WP)    :: Q
      real(WP)    :: X
      real(WP)    :: W
      real(WP)    :: sigma
      complex(WP) :: Delta
      real(WP)    :: E
      complex(WP) :: A
      complex(WP) :: dxi_r
      complex(WP) :: dxi_h
      complex(WP) :: deul_phi
      complex(WP) :: dF
      real(WP)    :: dtau

      ! Initialize the responses at this point

      rs%xi_r(j) = 0._WP
      rs%xi_h(j) = 0._WP
      rs%eul_phi(j) = 0._WP
      rs%F(j) = 0._WP
      rs%tau(j) = 0._WP

      ! Loop through the summary/detail data, adding in contributions

      n = dd(1)%n

      do i = 1, sd%n

         omega = sd%omega(i)%re
         gamma = -sd%omega(i)%im

         ! Calculate the orbital separation parameter

         R_a = (Omega_orb**2/(1._WP+or_p%q))**(1._WP/3._WP)

         ! Evaluate the various different terms involved in the
         ! Burkart+2012 formalism

         ! B+12, eqn. 8

         eps = or_p%q*(R_a/(1._WP-or_p%e))**(md_p%l+1)

         ! B+12, eqn. 9

         Q = -(2*md_p%l+1)/(4._WP*PI)*dd(i)%eul_phi(n)

         ! B+12, eqns. 11 & 12

         X = hansen_X(or_p%e, -(md_p%l+1), -md_p%m, -or_p%k)*(1._WP-or_p%e)**(md_p%l+1)

         ! B+12, eqn. A3

         W = 4._WP*PI/(2*md_p%l+1)*REAL(CONJG(spherical_Y(md_p%l, md_p%m, HALFPI, 0._WP)))

         ! B+12, eqn. 13

         sigma = or_p%k*Omega_orb - md_p%m*Omega_rot

         Delta = omega**2/((omega**2 - sigma**2) - 2._WP*CMPLX(0._WP, 1._WP, KIND=WP)*gamma*sigma)

         ! B+12, eqn.

         E = omega**2/(2._WP*PI)*dd(i)%E

         ! B+12, eqn. 7

         A = 2._WP*eps*Q*X*W*Delta/E

         ! Evaluate contributions

         dxi_r = A*dd(i)%xi_r(n)
         dxi_h = A*dd(i)%xi_h(n)
         deul_phi = A*dd(i)%eul_phi(n)

         dF = -Q*Delta*dd(i)%eul_phi(n)/E

         ! B+12, eqn C4 (with some transformation)

         dtau = 8._WP*(or_p%q*R_a**(md_p%l+1)*W*X*Q*ABS(Delta))**2 * md_p%m*sigma*gamma/E
         !dtau = 8._WP*or_p%q*R_a**(2*md_p%l+2)*W**2*X**2*Q**2 * md_p%m * omega**2*sigma*gamma/((omega**2 - sigma**2)**2 + 4*gamma**2*sigma**2) / E

         ! Add them in

         rs%xi_r(j) = rs%xi_r(j) + dxi_r
         rs%xi_h(j) = rs%xi_h(j) + dxi_h
         rs%eul_phi(j) = rs%eul_phi(j) + deul_phi
         rs%tau(j) = rs%tau(j) + dtau
         rs%F(j) = rs%F(j) + dF
         
         if (dd(i)%n_pg == rs%n_pg) then

            rs%xi_r_res(j) = dxi_r
            rs%xi_h_res(j) = dxi_h
            rs%eul_phi_res(j) = deul_phi
            rs%tau_res(j) = dtau
            rs%F_res(j) = dF

         endif

      end do

      ! Finish

      return

   end subroutine eval_response

   !****

   subroutine write_response(hg, rs)

      type(hgroup_t), intent(inout) :: hg
      type(response_t), intent(in)  :: rs

      ! Write the response_t

      call write_attr(hg, 'omega', rs%omega)
      call write_attr(hg, 'n_pg', rs%n_pg)

      call write_dset(hg, 'sigma', rs%sigma)

      call write_dset(hg, 'Omega_orb', rs%Omega_orb)
      call write_dset(hg, 'Omega_rot', rs%Omega_rot)

      call write_dset(hg, 'xi_r', rs%xi_r)
      call write_dset(hg, 'xi_h', rs%xi_h)
      call write_dset(hg, 'eul_phi', rs%eul_phi)
      call write_dset(hg, 'F', rs%F)
      call write_dset(hg, 'tau', rs%tau)
      
      call write_dset(hg, 'xi_r_res', rs%xi_r_res)
      call write_dset(hg, 'xi_h_res', rs%xi_h_res)
      call write_dset(hg, 'eul_phi_res', rs%eul_phi_res)
      call write_dset(hg, 'tau_res', rs%tau_res)
      
      ! Finish

      return

   end subroutine write_response

   !****

   function rs_name(i) result(name)

      integer, intent(in) :: i
      character(256)      :: name

      ! Construct the group name

      write(name, 100) i
100   format('rs(',I0,')')

      ! Finish

   end function rs_name

end program
