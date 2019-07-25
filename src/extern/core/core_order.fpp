! Module   : core_order
! Purpose  : sorting/searching

$include 'core.inc'

module core_order

  ! Uses

  use core_kinds

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Interfaces

  interface sort_indices
     module procedure sort_indices_i_i4_
     module procedure sort_indices_i_i8_
     module procedure sort_indices_r_sp_
     module procedure sort_indices_r_dp_
     module procedure sort_indices_a_
  end interface sort_indices

  interface unique_indices
     module procedure unique_indices_i_i4_
     module procedure unique_indices_i_i8_
     module procedure unique_indices_r_sp_
     module procedure unique_indices_r_dp_
     module procedure unique_indices_a_
  end interface unique_indices

  interface locate
     module procedure locate_i_i4_
     module procedure locate_i_i8_
     module procedure locate_r_sp_
     module procedure locate_r_dp_
     module procedure locate_a_
     module procedure locate_uniform_r_sp_
     module procedure locate_uniform_r_dp_
  end interface locate

  interface index_1d
     module procedure index_1d_i_i4_
     module procedure index_1d_i_i8_
  end interface index_1d

  interface index_nd
     module procedure index_nd_i_i4_
     module procedure index_nd_i_i8_
  end interface index_nd

  ! Access specifiers

  private

  public :: sort_indices
  public :: unique_indices
  public :: index_1d
  public :: index_nd
  public :: locate

  ! Procedures

contains

  $define $SORT_INDICES $sub

  $local $INFIX $1
  $local $TYPE $2

  function sort_indices_${INFIX}_ (x, descend) result (indices)

    $TYPE, intent(in)             :: x(:)
    logical, optional, intent(in) :: descend
    integer                       :: indices(SIZE(x))

    integer :: N
    integer :: L(0:SIZE(x)+1)
    integer :: i
    integer :: p
    integer :: q
    integer :: s
    integer :: t

    ! Calculate the indices to sort x in increasing order using the
    ! List Merge Sort algorithm described by Knuth (1998, The Art of
    ! Computer Programming Vol. 3, Addison- Wesley, Massachusetts, 3rd
    ! edn.). Labels and variable names deliberately follow those of
    ! Knuth (apart from K() -> x())

    N = SIZE(x)

    ! Special casing if N is 0 or 1; Knuth assumes N >= 2

    if(N == 0) then
       return
    elseif(N == 1) then
       indices(1) = 1
       return
    end if

    ! L1 : Prepare two lists

    L(0) = 1
    L(1:N-2) = [(-(i+2), i=1,N-2)]
    L(N-1:N) = 0
    L(N+1) = 2

    ! L2 : Begin new pass

    L2_loop : do

       s = 0
       t = N + 1

       p = L(s)
       q = L(t)

       if(q == 0) exit L2_loop

       ! L3 : Compare K_p:K_q

       L3_loop : do

          if(x(p) > x(q)) then

             ! L6 : Advance q

             L(s) = SIGN(q, L(s))
             s = q
             q = L(q)

             if(q > 0) cycle L3_loop

             ! L7 : Complete the sublist

             L(s) = p
             s = t

             L7_loop : do

                t = p
                p = L(p)

                if(p <= 0) exit L7_loop

             end do L7_loop

          else

             ! L4 : Advance q

             L(s) = SIGN(p, L(s))
             s = p
             p = L(p)

             if(p > 0) cycle L3_loop

             ! L5 : Complete the sublist

             L(s) = q
             s = t

             L5_loop : do

                t = q
                q = L(q)

                if(q <= 0) exit L5_loop

             end do L5_loop

          end if

          ! L8 : End of pass?
          
          p = -p
          q = -q

          if(q == 0) then
             L(s) = SIGN(p, L(s))
             L(t) = 0
             exit L3_loop
          endif

       end do L3_loop

    end do L2_loop

    ! Use the link list to set up the indices array
    
    indices(1) = L(0)

    do i = 2,N
       indices(i) = L(indices(i-1))
    enddo

    ! If necessary, reverse the indices

    if(PRESENT(descend)) then
       if(descend) indices = indices(N:1:-1)
    endif

    ! Finish

    return

  end function sort_indices_${INFIX}_

  $endsub

  $SORT_INDICES(i_i4,integer(I4))
  $SORT_INDICES(i_i8,integer(I8))
  $SORT_INDICES(r_sp,real(SP))
  $SORT_INDICES(r_dp,real(DP))
  $SORT_INDICES(a,character(*))

!****

  $define $UNIQUE_INDICES $sub

  $local $INFIX $1
  $local $TYPE $2

  function unique_indices_${INFIX}_ (x, descend) result (indices)

    $TYPE, intent(in)             :: x(:)
    logical, optional, intent(in) :: descend
    integer, allocatable          :: indices(:)

    integer :: indices_(SIZE(x))
    integer :: n
    integer :: i

    ! Calculate the indices to sort x in increasing order and then
    ! select unique values

    indices_ = sort_indices_${INFIX}_(x)

    n = 1

    do i = 2,SIZE(x)
       if(x(indices_(i)) /= x(indices_(n))) then
          n = n + 1
          indices_(n) = indices_(i)
       endif
    end do

    indices = indices_(:n)

    $if($DEBUG)
    do i = 1,SIZE(indices)
       $ASSERT(COUNT(x(indices(i)) == x(indices)) == 1)
    end do
    $endif

    ! If necessary, reverse the indices

    if(PRESENT(descend)) then
       if(descend) indices = indices(n:1:-1)
    endif

    ! Finish

    return

  end function unique_indices_${INFIX}_

  $endsub

  $UNIQUE_INDICES(i_i4,integer(I4))
  $UNIQUE_INDICES(i_i8,integer(I8))
  $UNIQUE_INDICES(r_sp,real(SP))
  $UNIQUE_INDICES(r_dp,real(DP))
  $UNIQUE_INDICES(a,character(*))

!****

  $define $INDEX_1d $sub

  $local $INFIX $1
  $local $TYPE $2

  function index_1d_${INFIX}_ (i, s) result (i1)

    $TYPE, intent(in)   :: i(:)
    integer, intent(in) :: s(:)
    $TYPE               :: i1

    integer :: r

    $CHECK_BOUNDS(SIZE(s),SIZE(i))

    ! Calculate a 1d index from an nd index

    i1 = 0

    do r = 1,SIZE(i)
       i1 = i1 + (i(r)-1)*PRODUCT(s(:r-1))
    end do

    i1 = i1 + 1

    ! Finish

    return

  end function index_1d_${INFIX}_

  $endsub

  $INDEX_1d(i_i4,integer(I4))
  $INDEX_1d(i_i8,integer(I8))

!****

  $define $INDEX_ND $sub

  $local $INFIX $1
  $local $TYPE $2

  function index_nd_${INFIX}_ (i1, s) result (i)

    $TYPE, intent(in)   :: i1
    integer, intent(in) :: s(:)
    $TYPE               :: i(SIZE(s))

    integer :: r

    $CHECK_BOUNDS(SIZE(s),SIZE(i))

    ! Calculate an nd index from a 1d index

    i(1) = i1 - 1

    do r = SIZE(s),2,-1
       i(r) = i(1)/PRODUCT(s(:r-1))
       i(1) = i(1) - i(r)*PRODUCT(s(:r-1))
    end do

    i = i + 1

    ! Finish

    return

  end function index_nd_${INFIX}_

  $endsub

  $INDEX_ND(i_i4,integer(I4))
  $INDEX_ND(i_i8,integer(I8))

!****

  $define $LOCATE $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine locate_${INFIX}_ (x, x_loc, i_loc)

    $TYPE, intent(in)    :: x(:)
    $TYPE, intent(in)    :: x_loc
    integer, intent(out) :: i_loc

    integer       :: n
    integer, save :: i_lo
    integer, save :: i_hi
    integer       :: di
    integer       :: i_mid

    !$OMP THREADPRIVATE (i_lo, i_hi)

    ! Use a binary search to find where x_loc falls in x (assumed to
    ! be in ascending order); x(i_loc) <= x_loc <= x(i_loc+1)

    n = SIZE(x)

    if(x_loc == x(n)) then

       i_loc = n - 1

    elseif(x_loc == x(1)) then

       i_loc = 1

    else

       if(i_lo < 1 .OR. i_lo >= n) then

          i_lo = 0
          i_hi = n+1

       else

          di = 1

          if(x_loc >= x(i_lo)) then

             search_up_loop : do

                i_hi = i_lo + di

                if(i_hi > n) then
                   i_hi = n + 1
                   exit search_up_loop
                endif

                if(x_loc < x(i_hi)) exit search_up_loop

                i_lo = i_hi
                di = 2*di

             end do search_up_loop

          else

             search_down_loop : do

                i_hi = i_lo
                i_lo = i_hi - di

                if(i_lo < 1) then
                   i_lo = 0
                   exit search_down_loop
                endif

                if(x_loc >= x(i_lo)) exit search_down_loop

             end do search_down_loop

          endif

       endif

       refine_loop : do

          if(i_hi-i_lo <= 1) exit refine_loop

          i_mid = (i_hi + i_lo)/2
       
          if(x_loc >= x(i_mid)) then
             i_lo = i_mid
          else
             i_hi = i_mid
          endif

       end do refine_loop

       i_loc = i_lo

    endif

    ! Finish

    return

  end subroutine locate_${INFIX}_

  $endsub

  $LOCATE(i_i4,integer(I4))
  $LOCATE(i_i8,integer(I8))
  $LOCATE(r_sp,real(SP))
  $LOCATE(r_dp,real(DP))
  $LOCATE(a,character(*))

!****

  $define $LOCATE_UNIFORM $sub

  $local $INFIX $1
  $local $TYPE $2

  subroutine locate_uniform_${INFIX}_ (x_0, dx, x_loc, i_loc)

    $TYPE, intent(in)    :: x_0
    $TYPE, intent(in)    :: dx
    $TYPE, intent(in)    :: x_loc
    integer, intent(out) :: i_loc

    ! Determine i_loc so that x_0 + (i_loc-1)*dx <= x_loc < x_0 + i_loc*dx

    i_loc = FLOOR((x_loc - x_0)/dx) + 1

    if(x_0 + (i_loc-1)*dx > x_loc) i_loc = i_loc - 1
    if(x_0 + i_loc*dx <= x_loc) i_loc = i_loc + 1

    ! Finish

    return

  end subroutine locate_uniform_${INFIX}_

  $endsub

  $LOCATE_UNIFORM(r_sp,real(SP))
  $LOCATE_UNIFORM(r_dp,real(DP))

end module core_order
