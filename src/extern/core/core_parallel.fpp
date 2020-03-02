! Module   : core_parallel
! Purpose  : parallel support

$include 'core.inc'

module core_parallel

  ! Uses

  use core_kinds
  use core_order

  use ISO_FORTRAN_ENV

  $if($MPI)
  use MPI
  $endif

  $if($OMP)
  use omp_lib
  $endif

  ! No implicit typing

  implicit none

  ! Module variables

  integer, save, protected :: MPI_SIZE
  integer, save, protected :: MPI_RANK
  integer, save, protected :: OMP_SIZE_MAX

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_i_i4_0_
     module procedure bcast_i_i4_1_
     module procedure bcast_i_i4_2_
     module procedure bcast_i_i4_3_
     module procedure bcast_i_i4_4_
     module procedure bcast_i_i8_0_
     module procedure bcast_i_i8_1_
     module procedure bcast_i_i8_2_
     module procedure bcast_i_i8_3_
     module procedure bcast_i_i8_4_
     module procedure bcast_r_sp_0_
     module procedure bcast_r_sp_1_
     module procedure bcast_r_sp_2_    
     module procedure bcast_r_sp_3_     
     module procedure bcast_r_sp_4_
     module procedure bcast_r_dp_0_
     module procedure bcast_r_dp_1_
     module procedure bcast_r_dp_2_     
     module procedure bcast_r_dp_3_
     module procedure bcast_r_dp_4_
     module procedure bcast_c_sp_0_
     module procedure bcast_c_sp_1_
     module procedure bcast_c_sp_2_   
     module procedure bcast_c_sp_3_    
     module procedure bcast_c_sp_4_
     module procedure bcast_c_dp_0_
     module procedure bcast_c_dp_1_
     module procedure bcast_c_dp_2_
     module procedure bcast_c_dp_3_
     module procedure bcast_c_dp_4_
     module procedure bcast_l_0_
     module procedure bcast_l_1_
     module procedure bcast_l_2_
     module procedure bcast_l_3_
     module procedure bcast_l_4_
     module procedure bcast_a_0_
     module procedure bcast_a_1_
     module procedure bcast_a_2_
     module procedure bcast_a_3_
     module procedure bcast_a_4_
  end interface bcast

  interface bcast_seq
     module procedure bcast_seq_i_i4_1_
     module procedure bcast_seq_i_i4_2_
     module procedure bcast_seq_i_i4_3_
     module procedure bcast_seq_i_i4_4_
     module procedure bcast_seq_i_i8_1_
     module procedure bcast_seq_i_i8_2_
     module procedure bcast_seq_i_i8_3_
     module procedure bcast_seq_i_i8_4_
     module procedure bcast_seq_r_sp_1_
     module procedure bcast_seq_r_sp_2_
     module procedure bcast_seq_r_sp_3_    
     module procedure bcast_seq_r_sp_4_
     module procedure bcast_seq_r_dp_1_
     module procedure bcast_seq_r_dp_2_
     module procedure bcast_seq_r_dp_3_
     module procedure bcast_seq_r_dp_4_
     module procedure bcast_seq_c_sp_1_
     module procedure bcast_seq_c_sp_2_
     module procedure bcast_seq_c_sp_3_    
     module procedure bcast_seq_c_sp_4_
     module procedure bcast_seq_c_dp_1_
     module procedure bcast_seq_c_dp_2_    
     module procedure bcast_seq_c_dp_3_    
     module procedure bcast_seq_c_dp_4_
     module procedure bcast_seq_l_1_
     module procedure bcast_seq_l_2_
     module procedure bcast_seq_l_3_
     module procedure bcast_seq_l_4_
  end interface bcast_seq

  interface bcast_alloc
     module procedure bcast_alloc_i_i4_0_
     module procedure bcast_alloc_i_i4_1_
     module procedure bcast_alloc_i_i4_2_
     module procedure bcast_alloc_i_i4_3_
     module procedure bcast_alloc_i_i4_4_
     module procedure bcast_alloc_i_i8_0_
     module procedure bcast_alloc_i_i8_1_
     module procedure bcast_alloc_i_i8_2_
     module procedure bcast_alloc_i_i8_3_
     module procedure bcast_alloc_i_i8_4_
     module procedure bcast_alloc_r_sp_0_
     module procedure bcast_alloc_r_sp_1_
     module procedure bcast_alloc_r_sp_2_     
     module procedure bcast_alloc_r_sp_3_
     module procedure bcast_alloc_r_sp_4_
     module procedure bcast_alloc_r_dp_0_
     module procedure bcast_alloc_r_dp_1_
     module procedure bcast_alloc_r_dp_2_
     module procedure bcast_alloc_r_dp_3_
     module procedure bcast_alloc_r_dp_4_
     module procedure bcast_alloc_c_sp_0_
     module procedure bcast_alloc_c_sp_1_
     module procedure bcast_alloc_c_sp_2_
     module procedure bcast_alloc_c_sp_3_
     module procedure bcast_alloc_c_sp_4_
     module procedure bcast_alloc_c_dp_0_
     module procedure bcast_alloc_c_dp_1_
     module procedure bcast_alloc_c_dp_2_
     module procedure bcast_alloc_c_dp_3_
     module procedure bcast_alloc_c_dp_4_
     module procedure bcast_alloc_l_0_
     module procedure bcast_alloc_l_1_
     module procedure bcast_alloc_l_2_
     module procedure bcast_alloc_l_3_
     module procedure bcast_alloc_l_4_
     module procedure bcast_alloc_a_0_
     module procedure bcast_alloc_a_1_
     module procedure bcast_alloc_a_2_
     module procedure bcast_alloc_a_3_
     module procedure bcast_alloc_a_4_
  end interface bcast_alloc

  interface send
     module procedure send_i_i4_0_
     module procedure send_i_i4_1_
     module procedure send_i_i4_2_
     module procedure send_i_i4_3_
     module procedure send_i_i4_4_
     module procedure send_i_i8_0_
     module procedure send_i_i8_1_
     module procedure send_i_i8_2_
     module procedure send_i_i8_3_
     module procedure send_i_i8_4_
     module procedure send_r_sp_0_
     module procedure send_r_sp_1_
     module procedure send_r_sp_2_     
     module procedure send_r_sp_3_     
     module procedure send_r_sp_4_
     module procedure send_r_dp_0_
     module procedure send_r_dp_1_
     module procedure send_r_dp_2_     
     module procedure send_r_dp_3_     
     module procedure send_r_dp_4_
     module procedure send_c_sp_0_
     module procedure send_c_sp_1_
     module procedure send_c_sp_2_     
     module procedure send_c_sp_3_     
     module procedure send_c_sp_4_
     module procedure send_c_dp_0_
     module procedure send_c_dp_1_
     module procedure send_c_dp_2_     
     module procedure send_c_dp_3_     
     module procedure send_c_dp_4_
     module procedure send_l_0_
     module procedure send_l_1_
     module procedure send_l_2_
     module procedure send_l_3_
     module procedure send_l_4_
  end interface send

  interface recv
     module procedure recv_i_i4_0_
     module procedure recv_i_i4_1_
     module procedure recv_i_i4_2_
     module procedure recv_i_i4_3_
     module procedure recv_i_i4_4_
     module procedure recv_i_i8_0_
     module procedure recv_i_i8_1_
     module procedure recv_i_i8_2_
     module procedure recv_i_i8_3_
     module procedure recv_i_i8_4_
     module procedure recv_r_sp_0_
     module procedure recv_r_sp_1_
     module procedure recv_r_sp_2_     
     module procedure recv_r_sp_3_     
     module procedure recv_r_sp_4_
     module procedure recv_r_dp_0_
     module procedure recv_r_dp_1_
     module procedure recv_r_dp_2_     
     module procedure recv_r_dp_3_     
     module procedure recv_r_dp_4_
     module procedure recv_c_sp_0_
     module procedure recv_c_sp_1_
     module procedure recv_c_sp_2_     
     module procedure recv_c_sp_3_     
     module procedure recv_c_sp_4_
     module procedure recv_c_dp_0_
     module procedure recv_c_dp_1_
     module procedure recv_c_dp_2_     
     module procedure recv_c_dp_3_     
     module procedure recv_c_dp_4_
     module procedure recv_l_0_
     module procedure recv_l_1_
     module procedure recv_l_2_
     module procedure recv_l_3_
     module procedure recv_l_4_
  end interface recv

  interface recv_any
     module procedure recv_any_i_i4_0_
     module procedure recv_any_i_i4_1_
     module procedure recv_any_i_i4_2_
     module procedure recv_any_i_i4_3_
     module procedure recv_any_i_i4_4_
     module procedure recv_any_i_i8_0_
     module procedure recv_any_i_i8_1_
     module procedure recv_any_i_i8_2_
     module procedure recv_any_i_i8_3_
     module procedure recv_any_i_i8_4_
     module procedure recv_any_r_sp_0_
     module procedure recv_any_r_sp_1_
     module procedure recv_any_r_sp_2_     
     module procedure recv_any_r_sp_3_     
     module procedure recv_any_r_sp_4_
     module procedure recv_any_r_dp_0_
     module procedure recv_any_r_dp_1_
     module procedure recv_any_r_dp_2_     
     module procedure recv_any_r_dp_3_     
     module procedure recv_any_r_dp_4_
     module procedure recv_any_c_sp_0_
     module procedure recv_any_c_sp_1_
     module procedure recv_any_c_sp_2_     
     module procedure recv_any_c_sp_3_     
     module procedure recv_any_c_sp_4_
     module procedure recv_any_c_dp_0_
     module procedure recv_any_c_dp_1_
     module procedure recv_any_c_dp_2_     
     module procedure recv_any_c_dp_3_     
     module procedure recv_any_c_dp_4_
     module procedure recv_any_l_0_
     module procedure recv_any_l_1_
     module procedure recv_any_l_2_
     module procedure recv_any_l_3_
     module procedure recv_any_l_4_
  end interface recv_any

  interface gatherv
     module procedure gatherv_i_i4_0_
     module procedure gatherv_i_i4_1_
     module procedure gatherv_i_i4_2_
     module procedure gatherv_i_i4_3_
     module procedure gatherv_i_i4_4_
     module procedure gatherv_i_i8_0_
     module procedure gatherv_i_i8_1_
     module procedure gatherv_i_i8_2_
     module procedure gatherv_i_i8_3_
     module procedure gatherv_i_i8_4_
     module procedure gatherv_r_sp_0_
     module procedure gatherv_r_sp_1_
     module procedure gatherv_r_sp_2_     
     module procedure gatherv_r_sp_3_     
     module procedure gatherv_r_sp_4_
     module procedure gatherv_r_dp_0_
     module procedure gatherv_r_dp_1_
     module procedure gatherv_r_dp_2_     
     module procedure gatherv_r_dp_3_     
     module procedure gatherv_r_dp_4_
     module procedure gatherv_c_sp_0_
     module procedure gatherv_c_sp_1_
     module procedure gatherv_c_sp_2_     
     module procedure gatherv_c_sp_3_     
     module procedure gatherv_c_sp_4_
     module procedure gatherv_c_dp_0_
     module procedure gatherv_c_dp_1_
     module procedure gatherv_c_dp_2_     
     module procedure gatherv_c_dp_3_     
     module procedure gatherv_c_dp_4_
  end interface gatherv

  interface allgatherv
     module procedure allgatherv_i_i4_0_
     module procedure allgatherv_i_i4_1_
     module procedure allgatherv_i_i4_2_
     module procedure allgatherv_i_i4_3_
     module procedure allgatherv_i_i4_4_
     module procedure allgatherv_i_i8_0_
     module procedure allgatherv_i_i8_1_
     module procedure allgatherv_i_i8_2_
     module procedure allgatherv_i_i8_3_
     module procedure allgatherv_i_i8_4_
     module procedure allgatherv_r_sp_0_
     module procedure allgatherv_r_sp_1_
     module procedure allgatherv_r_sp_2_     
     module procedure allgatherv_r_sp_3_     
     module procedure allgatherv_r_sp_4_
     module procedure allgatherv_r_dp_0_
     module procedure allgatherv_r_dp_1_
     module procedure allgatherv_r_dp_2_     
     module procedure allgatherv_r_dp_3_     
     module procedure allgatherv_r_dp_4_
     module procedure allgatherv_c_sp_0_
     module procedure allgatherv_c_sp_1_
     module procedure allgatherv_c_sp_2_     
     module procedure allgatherv_c_sp_3_     
     module procedure allgatherv_c_sp_4_
     module procedure allgatherv_c_dp_0_
     module procedure allgatherv_c_dp_1_
     module procedure allgatherv_c_dp_2_     
     module procedure allgatherv_c_dp_3_     
     module procedure allgatherv_c_dp_4_
  end interface allgatherv

  interface allreduce
     module procedure allreduce_i_i4_0_
     module procedure allreduce_i_i4_1_
     module procedure allreduce_i_i4_2_
     module procedure allreduce_i_i4_3_
     module procedure allreduce_i_i4_4_
     module procedure allreduce_i_i8_0_
     module procedure allreduce_i_i8_1_
     module procedure allreduce_i_i8_2_
     module procedure allreduce_i_i8_3_
     module procedure allreduce_i_i8_4_
     module procedure allreduce_r_sp_0_
     module procedure allreduce_r_sp_1_
     module procedure allreduce_r_sp_2_     
     module procedure allreduce_r_sp_3_     
     module procedure allreduce_r_sp_4_
     module procedure allreduce_r_dp_0_
     module procedure allreduce_r_dp_1_
     module procedure allreduce_r_dp_2_     
     module procedure allreduce_r_dp_3_     
     module procedure allreduce_r_dp_4_
     module procedure allreduce_c_sp_0_
     module procedure allreduce_c_sp_1_
     module procedure allreduce_c_sp_2_     
     module procedure allreduce_c_sp_3_     
     module procedure allreduce_c_sp_4_
     module procedure allreduce_c_dp_0_
     module procedure allreduce_c_dp_1_
     module procedure allreduce_c_dp_2_     
     module procedure allreduce_c_dp_3_     
     module procedure allreduce_c_dp_4_
  end interface allreduce

  $endif

  ! Access specifiers

  private

  public :: MPI_SIZE
  public :: MPI_RANK
  public :: OMP_SIZE_MAX
  public :: init_parallel
  public :: final_parallel
  public :: omp_size
  public :: omp_rank
  $if($MPI)
  public :: MPI_COMM_WORLD
  public :: MPI_SUM
  public :: barrier
  public :: bcast
  public :: bcast_seq
  public :: bcast_alloc
  public :: send
  public :: recv
  public :: recv_any
  public :: gatherv
  public :: allgatherv
  public :: allreduce
  $endif
  $if($OMP)
  public :: omp_get_thread_num
  $endif
  public :: partition_tasks

contains

  subroutine init_parallel ()

    $if($MPI)
    integer :: mpi_err
    $endif

    ! Initialize MPI

    $if($MPI)

    call MPI_INIT(mpi_err)
    $ASSERT(mpi_err == MPI_SUCCESS,MPI initialization failed)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_SIZE, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, mpi_err)

    $else

    MPI_SIZE = 1
    MPI_RANK = 0

    $endif

    ! Initialize OpenMP

    $if($OMP)
    OMP_SIZE_MAX = omp_get_max_threads()
    $else
    OMP_SIZE_MAX = 1
    $endif

    ! Finish

    return

  end subroutine init_parallel

!****

  subroutine final_parallel ()
 
    $if($MPI)
    integer :: mpi_err
    $endif

    ! Finalize MPI

    $if($MPI)
    
    call MPI_FINALIZE(mpi_err)
    $ASSERT(mpi_err == MPI_SUCCESS,MPI finalization failed)

    $endif

    MPI_SIZE = 0
    MPI_RANK = 0

    ! Finish

    return

  end subroutine final_parallel

!****

  function omp_size ()

    integer :: omp_size

    ! Get the OpenMP thread size

    $if($OMP)
    omp_size = omp_get_num_threads()
    $else
    omp_size = 1
    $endif

    ! Finish

    return

  end function omp_size

!****

  function omp_rank ()

    integer :: omp_rank

    ! Get the OpenMP thread rank

    $if($OMP)
    omp_rank = omp_get_thread_num()
    $else
    omp_rank = 0
    $endif

    ! finish

    return

  end function omp_rank

!****

  $if($MPI)

  subroutine barrier ()

    integer :: mpi_err

    ! Set up a barrier

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    $ASSERT(mpi_err == MPI_SUCCESS,MPI_BARRIER failed)

    ! Finish

    return

  end subroutine barrier

!****

  $define $BCAST $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_MPI_TYPE $3
  $local $BUFFER_RANK $4

  subroutine bcast_${INFIX}_${BUFFER_RANK}_ (buffer, root_rank)

    $BUFFER_TYPE, intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: root_rank

    integer :: mpi_err

    ! Broadcast the buffer

    call MPI_BCAST(buffer, PRODUCT(SHAPE(buffer)), $BUFFER_MPI_TYPE, root_rank, MPI_COMM_WORLD, mpi_err)
    $ASSERT(mpi_err == MPI_SUCCESS,MPI_BCAST failed)
    
    ! This barrier is required due to an issue with MVAPICH
    ! (non-master ranks do not return from BCAST until another MPI
    ! operation is issued)

    call barrier()

    ! Finish

    return

  end subroutine bcast_${INFIX}_${BUFFER_RANK}_

  $endsub

  $BCAST(i_i4,integer(I4),MPI_INTEGER4,0)
  $BCAST(i_i4,integer(I4),MPI_INTEGER4,1)
  $BCAST(i_i4,integer(I4),MPI_INTEGER4,2)
  $BCAST(i_i4,integer(I4),MPI_INTEGER4,3)
  $BCAST(i_i4,integer(I4),MPI_INTEGER4,4)

  $BCAST(i_i8,integer(I8),MPI_INTEGER8,0)
  $BCAST(i_i8,integer(I8),MPI_INTEGER8,1)
  $BCAST(i_i8,integer(I8),MPI_INTEGER8,2)
  $BCAST(i_i8,integer(I8),MPI_INTEGER8,3)
  $BCAST(i_i8,integer(I8),MPI_INTEGER8,4)

  $BCAST(r_sp,real(SP),MPI_REAL,0)
  $BCAST(r_sp,real(SP),MPI_REAL,1)
  $BCAST(r_sp,real(SP),MPI_REAL,2)
  $BCAST(r_sp,real(SP),MPI_REAL,3)
  $BCAST(r_sp,real(SP),MPI_REAL,4)

  $BCAST(r_dp,real(DP),MPI_DOUBLE_PRECISION,0)
  $BCAST(r_dp,real(DP),MPI_DOUBLE_PRECISION,1)
  $BCAST(r_dp,real(DP),MPI_DOUBLE_PRECISION,2)
  $BCAST(r_dp,real(DP),MPI_DOUBLE_PRECISION,3)
  $BCAST(r_dp,real(DP),MPI_DOUBLE_PRECISION,4)

  $BCAST(c_sp,complex(SP),MPI_COMPLEX,0)
  $BCAST(c_sp,complex(SP),MPI_COMPLEX,1)
  $BCAST(c_sp,complex(SP),MPI_COMPLEX,2)
  $BCAST(c_sp,complex(SP),MPI_COMPLEX,3)
  $BCAST(c_sp,complex(SP),MPI_COMPLEX,4)

  $BCAST(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,0)
  $BCAST(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,1)
  $BCAST(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,2)
  $BCAST(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,3)
  $BCAST(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,4)

  $BCAST(l,logical,MPI_LOGICAL,0)
  $BCAST(l,logical,MPI_LOGICAL,1)
  $BCAST(l,logical,MPI_LOGICAL,2)
  $BCAST(l,logical,MPI_LOGICAL,3)
  $BCAST(l,logical,MPI_LOGICAL,4)

!****

  $define $BCAST_SPECIAL $sub

  $local $BUFFER_RANK $1

  subroutine bcast_a_${BUFFER_RANK}_ (buffer, root)

    character(*), intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: root

    integer :: mpi_err

    ! Broadcast the buffer

    call MPI_BCAST(buffer, LEN(buffer)*PRODUCT(SHAPE(buffer)), MPI_CHARACTER, root, MPI_COMM_WORLD, mpi_err)
    $ASSERT(mpi_err == MPI_SUCCESS,MPI_BCAST failed)

    ! This barrier is required due to an issue with MVAPICH
    ! (non-master ranks do not return from BCAST until another MPI
    ! operation is issued)

    call barrier()

    ! Finish

    return

  end subroutine bcast_a_${BUFFER_RANK}_

  $endsub

  $BCAST_SPECIAL(0)
  $BCAST_SPECIAL(1)
  $BCAST_SPECIAL(2)
  $BCAST_SPECIAL(3)
  $BCAST_SPECIAL(4)

!****

  $define $BCAST_SEQ $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_MPI_TYPE $3
  $local $BUFFER_RANK $4

  subroutine bcast_seq_${INFIX}_${BUFFER_RANK}_ (buffer, i1_a, i1_b, root_rank)

    $BUFFER_TYPE, intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: i1_b
    integer, intent(in)         :: i1_a
    integer, intent(in)         :: root_rank

    integer :: i_a($BUFFER_RANK)
    integer :: n
    integer :: mpi_err

    ! Broadcast the buffer sequence extending from 1-D index i1_a to
    ! 1-D index i1_b

    i_a = index_nd(i1_a, SHAPE(buffer))

    n = i1_b - i1_a + 1

    call MPI_BCAST(buffer($ARRAY_EXPAND(i_a,$BUFFER_RANK)), n, $BUFFER_MPI_TYPE, root_rank, MPI_COMM_WORLD, mpi_err)
    $ASSERT(mpi_err == MPI_SUCCESS,MPI_BCAST failed)
    
    ! This barrier is required due to an issue with MVAPICH
    ! (non-master ranks do not return from BCAST until another MPI
    ! operation is issued)

    call barrier()

    ! Finish

    return

  end subroutine bcast_seq_${INFIX}_${BUFFER_RANK}_

  $endsub

  $BCAST_SEQ(i_i4,integer(I4),MPI_INTEGER4,1)
  $BCAST_SEQ(i_i4,integer(I4),MPI_INTEGER4,2)
  $BCAST_SEQ(i_i4,integer(I4),MPI_INTEGER4,3)
  $BCAST_SEQ(i_i4,integer(I4),MPI_INTEGER4,4)

  $BCAST_SEQ(i_i8,integer(I8),MPI_INTEGER8,1)
  $BCAST_SEQ(i_i8,integer(I8),MPI_INTEGER8,2)
  $BCAST_SEQ(i_i8,integer(I8),MPI_INTEGER8,3)
  $BCAST_SEQ(i_i8,integer(I8),MPI_INTEGER8,4)

  $BCAST_SEQ(r_sp,real(SP),MPI_REAL,1)
  $BCAST_SEQ(r_sp,real(SP),MPI_REAL,2)
  $BCAST_SEQ(r_sp,real(SP),MPI_REAL,3)
  $BCAST_SEQ(r_sp,real(SP),MPI_REAL,4)

  $BCAST_SEQ(r_dp,real(DP),MPI_DOUBLE_PRECISION,1)
  $BCAST_SEQ(r_dp,real(DP),MPI_DOUBLE_PRECISION,2)
  $BCAST_SEQ(r_dp,real(DP),MPI_DOUBLE_PRECISION,3)
  $BCAST_SEQ(r_dp,real(DP),MPI_DOUBLE_PRECISION,4)

  $BCAST_SEQ(c_sp,complex(SP),MPI_COMPLEX,1)
  $BCAST_SEQ(c_sp,complex(SP),MPI_COMPLEX,2)
  $BCAST_SEQ(c_sp,complex(SP),MPI_COMPLEX,3)
  $BCAST_SEQ(c_sp,complex(SP),MPI_COMPLEX,4)

  $BCAST_SEQ(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,1)
  $BCAST_SEQ(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,2)
  $BCAST_SEQ(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,3)
  $BCAST_SEQ(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,4)

  $BCAST_SEQ(l,logical,MPI_LOGICAL,1)
  $BCAST_SEQ(l,logical,MPI_LOGICAL,2)
  $BCAST_SEQ(l,logical,MPI_LOGICAL,3)
  $BCAST_SEQ(l,logical,MPI_LOGICAL,4)

!****

  $define $BCAST_ALLOC $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_RANK $3

  subroutine bcast_alloc_${INFIX}_${BUFFER_RANK}_ (buffer, root_rank)

    $BUFFER_TYPE, allocatable, intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)                      :: root_rank

    logical :: alloc
    $if ($BUFFER_RANK > 0)
    integer :: lb($BUFFER_RANK)
    integer :: ub($BUFFER_RANK)
    $endif

    ! Deallocate the buffer on non-root processors

    if (MPI_RANK /= root_rank .AND. ALLOCATED(buffer)) then
       deallocate(buffer)
    endif

    ! Check if the buffer is allocated on the root processor
    
    if (MPI_RANK == root_rank) alloc = ALLOCATED(buffer)
    call bcast(alloc, root_rank)

    if (alloc) then

       ! Broadcast the buffer bounds

       $if ($BUFFER_RANK > 0)

       if (MPI_RANK == root_rank) then
          lb = LBOUND(buffer)
          ub = UBOUND(buffer)
       endif

       call bcast(lb, root_rank)
       call bcast(ub, root_rank)

       $endif

       ! Allocate the buffer

       $if ($BUFFER_RANK > 0)

       if (MPI_RANK /= root_rank) allocate(buffer($ARRAY_EXPAND(lb,ub,$BUFFER_RANK)))

       $else

       if (MPI_RANK /= root_rank) allocate(buffer)

       $endif

       ! Broadcast the buffer

       call bcast(buffer, root_rank)

    endif

    ! Finish

    return

  end subroutine bcast_alloc_${INFIX}_${BUFFER_RANK}_

  $endsub

  $BCAST_ALLOC(i_i4,integer(I4),0)
  $BCAST_ALLOC(i_i4,integer(I4),1)
  $BCAST_ALLOC(i_i4,integer(I4),2)
  $BCAST_ALLOC(i_i4,integer(I4),3)
  $BCAST_ALLOC(i_i4,integer(I4),4)

  $BCAST_ALLOC(i_i8,integer(I8),0)
  $BCAST_ALLOC(i_i8,integer(I8),1)
  $BCAST_ALLOC(i_i8,integer(I8),2)
  $BCAST_ALLOC(i_i8,integer(I8),3)
  $BCAST_ALLOC(i_i8,integer(I8),4)

  $BCAST_ALLOC(r_sp,real(SP),0)
  $BCAST_ALLOC(r_sp,real(SP),1)
  $BCAST_ALLOC(r_sp,real(SP),2)
  $BCAST_ALLOC(r_sp,real(SP),3)
  $BCAST_ALLOC(r_sp,real(SP),4)

  $BCAST_ALLOC(r_dp,real(DP),0)
  $BCAST_ALLOC(r_dp,real(DP),1)
  $BCAST_ALLOC(r_dp,real(DP),2)
  $BCAST_ALLOC(r_dp,real(DP),3)
  $BCAST_ALLOC(r_dp,real(DP),4)

  $BCAST_ALLOC(c_sp,complex(SP),0)
  $BCAST_ALLOC(c_sp,complex(SP),1)
  $BCAST_ALLOC(c_sp,complex(SP),2)
  $BCAST_ALLOC(c_sp,complex(SP),3)
  $BCAST_ALLOC(c_sp,complex(SP),4)

  $BCAST_ALLOC(c_dp,complex(DP),0)
  $BCAST_ALLOC(c_dp,complex(DP),1)
  $BCAST_ALLOC(c_dp,complex(DP),2)
  $BCAST_ALLOC(c_dp,complex(DP),3)
  $BCAST_ALLOC(c_dp,complex(DP),4)

  $BCAST_ALLOC(a,character(*),0)
  $BCAST_ALLOC(a,character(*),1)
  $BCAST_ALLOC(a,character(*),2)
  $BCAST_ALLOC(a,character(*),3)
  $BCAST_ALLOC(a,character(*),4)

  $BCAST_ALLOC(l,logical,0)
  $BCAST_ALLOC(l,logical,1)
  $BCAST_ALLOC(l,logical,2)
  $BCAST_ALLOC(l,logical,3)
  $BCAST_ALLOC(l,logical,4)

!****

  $define $SEND $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_MPI_TYPE $3
  $local $BUFFER_RANK $4

  subroutine send_${INFIX}_${BUFFER_RANK}_ (buffer, dest_rank, tag, sync)

    $BUFFER_TYPE, intent(in)      :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)           :: dest_rank
    integer, intent(in), optional :: tag
    logical, intent(in), optional :: sync

    integer :: tag_
    logical :: sync_
    integer :: mpi_err

    if(PRESENT(tag)) then
       tag_ = tag
    else
       tag_ = 0
    endif

    if(PRESENT(sync)) then
       sync_ = sync
    else
       sync_ = .FALSE.
    endif

    ! Send the buffer

    if(sync_) then
       call MPI_SSEND(buffer, PRODUCT(SHAPE(buffer)), $BUFFER_MPI_TYPE, dest_rank, tag_, MPI_COMM_WORLD, mpi_err)
       $ASSERT(mpi_err == MPI_SUCCESS,MPI_SSEND failed)
    else
       call MPI_SEND(buffer, PRODUCT(SHAPE(buffer)), $BUFFER_MPI_TYPE, dest_rank, tag_, MPI_COMM_WORLD, mpi_err)
       $ASSERT(mpi_err == MPI_SUCCESS,MPI_SEND failed)
    endif
    
    ! Finish

    return

  end subroutine send_${INFIX}_${BUFFER_RANK}_

  $endsub

  $SEND(i_i4,integer(I4),MPI_INTEGER4,0)
  $SEND(i_i4,integer(I4),MPI_INTEGER4,1)
  $SEND(i_i4,integer(I4),MPI_INTEGER4,2)
  $SEND(i_i4,integer(I4),MPI_INTEGER4,3)
  $SEND(i_i4,integer(I4),MPI_INTEGER4,4)

  $SEND(i_i8,integer(I8),MPI_INTEGER8,0)
  $SEND(i_i8,integer(I8),MPI_INTEGER8,1)
  $SEND(i_i8,integer(I8),MPI_INTEGER8,2)
  $SEND(i_i8,integer(I8),MPI_INTEGER8,3)
  $SEND(i_i8,integer(I8),MPI_INTEGER8,4)

  $SEND(r_sp,real(SP),MPI_REAL,0)
  $SEND(r_sp,real(SP),MPI_REAL,1)
  $SEND(r_sp,real(SP),MPI_REAL,2)
  $SEND(r_sp,real(SP),MPI_REAL,3)
  $SEND(r_sp,real(SP),MPI_REAL,4)

  $SEND(r_dp,real(DP),MPI_DOUBLE_PRECISION,0)
  $SEND(r_dp,real(DP),MPI_DOUBLE_PRECISION,1)
  $SEND(r_dp,real(DP),MPI_DOUBLE_PRECISION,2)
  $SEND(r_dp,real(DP),MPI_DOUBLE_PRECISION,3)
  $SEND(r_dp,real(DP),MPI_DOUBLE_PRECISION,4)

  $SEND(c_sp,complex(SP),MPI_COMPLEX,0)
  $SEND(c_sp,complex(SP),MPI_COMPLEX,1)
  $SEND(c_sp,complex(SP),MPI_COMPLEX,2)
  $SEND(c_sp,complex(SP),MPI_COMPLEX,3)
  $SEND(c_sp,complex(SP),MPI_COMPLEX,4)

  $SEND(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,0)
  $SEND(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,1)
  $SEND(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,2)
  $SEND(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,3)
  $SEND(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,4)

  $SEND(l,logical,MPI_LOGICAL,0)
  $SEND(l,logical,MPI_LOGICAL,1)
  $SEND(l,logical,MPI_LOGICAL,2)
  $SEND(l,logical,MPI_LOGICAL,3)
  $SEND(l,logical,MPI_LOGICAL,4)

!****

  $define $RECV $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_MPI_TYPE $3
  $local $BUFFER_RANK $4

  subroutine recv_${INFIX}_${BUFFER_RANK}_ (buffer, src_rank, tag)

    $BUFFER_TYPE, intent(out)     :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)           :: src_rank
    integer, intent(in), optional :: tag

    integer :: tag_
    integer :: mpi_err

    if(PRESENT(tag)) then
       tag_ = tag
    else
       tag_ = 0
    endif

    ! Receive the buffer

    call MPI_RECV(buffer, PRODUCT(SHAPE(buffer)), $BUFFER_MPI_TYPE, src_rank, tag_, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    $ASSERT(mpi_err == MPI_SUCCESS,MPI_RECV failed)
    
    ! Finish

    return

  end subroutine recv_${INFIX}_${BUFFER_RANK}_

  $endsub

  $RECV(i_i4,integer(I4),MPI_INTEGER4,0)
  $RECV(i_i4,integer(I4),MPI_INTEGER4,1)
  $RECV(i_i4,integer(I4),MPI_INTEGER4,2)
  $RECV(i_i4,integer(I4),MPI_INTEGER4,3)
  $RECV(i_i4,integer(I4),MPI_INTEGER4,4)

  $RECV(i_i8,integer(I8),MPI_INTEGER8,0)
  $RECV(i_i8,integer(I8),MPI_INTEGER8,1)
  $RECV(i_i8,integer(I8),MPI_INTEGER8,2)
  $RECV(i_i8,integer(I8),MPI_INTEGER8,3)
  $RECV(i_i8,integer(I8),MPI_INTEGER8,4)

  $RECV(r_sp,real(SP),MPI_REAL,0)
  $RECV(r_sp,real(SP),MPI_REAL,1)
  $RECV(r_sp,real(SP),MPI_REAL,2)
  $RECV(r_sp,real(SP),MPI_REAL,3)
  $RECV(r_sp,real(SP),MPI_REAL,4)

  $RECV(r_dp,real(DP),MPI_DOUBLE_PRECISION,0)
  $RECV(r_dp,real(DP),MPI_DOUBLE_PRECISION,1)
  $RECV(r_dp,real(DP),MPI_DOUBLE_PRECISION,2)
  $RECV(r_dp,real(DP),MPI_DOUBLE_PRECISION,3)
  $RECV(r_dp,real(DP),MPI_DOUBLE_PRECISION,4)

  $RECV(c_sp,complex(SP),MPI_COMPLEX,0)
  $RECV(c_sp,complex(SP),MPI_COMPLEX,1)
  $RECV(c_sp,complex(SP),MPI_COMPLEX,2)
  $RECV(c_sp,complex(SP),MPI_COMPLEX,3)
  $RECV(c_sp,complex(SP),MPI_COMPLEX,4)

  $RECV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,0)
  $RECV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,1)
  $RECV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,2)
  $RECV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,3)
  $RECV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,4)

  $RECV(l,logical,MPI_LOGICAL,0)
  $RECV(l,logical,MPI_LOGICAL,1)
  $RECV(l,logical,MPI_LOGICAL,2)
  $RECV(l,logical,MPI_LOGICAL,3)
  $RECV(l,logical,MPI_LOGICAL,4)

!****

  $define $RECV_ANY $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_MPI_TYPE $3
  $local $BUFFER_RANK $4

  subroutine recv_any_${INFIX}_${BUFFER_RANK}_ (buffer, src_rank, tag)

    $BUFFER_TYPE, intent(out)     :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(out)          :: src_rank
    integer, intent(in), optional :: tag

    integer :: tag_
    integer :: status(MPI_STATUS_SIZE)
    integer :: mpi_err

    if(PRESENT(tag)) then
       tag_ = tag
    else
       tag_ = 0
    endif

    ! Receive the buffer

    call MPI_RECV(buffer, PRODUCT(SHAPE(buffer)), $BUFFER_MPI_TYPE, MPI_ANY_SOURCE, tag_, MPI_COMM_WORLD, status, mpi_err)
    $ASSERT(mpi_err == MPI_SUCCESS,MPI_RECV failed)

    src_rank = status(MPI_SOURCE)

    ! Finish

    return

  end subroutine recv_any_${INFIX}_${BUFFER_RANK}_

  $endsub

  $RECV_ANY(i_i4,integer(I4),MPI_INTEGER4,0)
  $RECV_ANY(i_i4,integer(I4),MPI_INTEGER4,1)
  $RECV_ANY(i_i4,integer(I4),MPI_INTEGER4,2)
  $RECV_ANY(i_i4,integer(I4),MPI_INTEGER4,3)
  $RECV_ANY(i_i4,integer(I4),MPI_INTEGER4,4)

  $RECV_ANY(i_i8,integer(I8),MPI_INTEGER8,0)
  $RECV_ANY(i_i8,integer(I8),MPI_INTEGER8,1)
  $RECV_ANY(i_i8,integer(I8),MPI_INTEGER8,2)
  $RECV_ANY(i_i8,integer(I8),MPI_INTEGER8,3)
  $RECV_ANY(i_i8,integer(I8),MPI_INTEGER8,4)

  $RECV_ANY(r_sp,real(SP),MPI_REAL,0)
  $RECV_ANY(r_sp,real(SP),MPI_REAL,1)
  $RECV_ANY(r_sp,real(SP),MPI_REAL,2)
  $RECV_ANY(r_sp,real(SP),MPI_REAL,3)
  $RECV_ANY(r_sp,real(SP),MPI_REAL,4)

  $RECV_ANY(r_dp,real(DP),MPI_DOUBLE_PRECISION,0)
  $RECV_ANY(r_dp,real(DP),MPI_DOUBLE_PRECISION,1)
  $RECV_ANY(r_dp,real(DP),MPI_DOUBLE_PRECISION,2)
  $RECV_ANY(r_dp,real(DP),MPI_DOUBLE_PRECISION,3)
  $RECV_ANY(r_dp,real(DP),MPI_DOUBLE_PRECISION,4)

  $RECV_ANY(c_sp,complex(SP),MPI_COMPLEX,0)
  $RECV_ANY(c_sp,complex(SP),MPI_COMPLEX,1)
  $RECV_ANY(c_sp,complex(SP),MPI_COMPLEX,2)
  $RECV_ANY(c_sp,complex(SP),MPI_COMPLEX,3)
  $RECV_ANY(c_sp,complex(SP),MPI_COMPLEX,4)

  $RECV_ANY(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,0)
  $RECV_ANY(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,1)
  $RECV_ANY(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,2)
  $RECV_ANY(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,3)
  $RECV_ANY(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,4)

  $RECV_ANY(l,logical,MPI_LOGICAL,0)
  $RECV_ANY(l,logical,MPI_LOGICAL,1)
  $RECV_ANY(l,logical,MPI_LOGICAL,2)
  $RECV_ANY(l,logical,MPI_LOGICAL,3)
  $RECV_ANY(l,logical,MPI_LOGICAL,4)

!****

  $define $GATHERV $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_MPI_TYPE $3
  $local $BUFFER_RANK $4

  subroutine gatherv_${INFIX}_${BUFFER_RANK}_ (send_buffer, sendcount, recv_buffer, recvcounts, displs, root_rank)

    $BUFFER_TYPE, intent(in)    :: send_buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: sendcount
    $BUFFER_TYPE, intent(inout) :: recv_buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: recvcounts(:)
    integer, intent(in)         :: displs(:)
    integer, intent(in)         :: root_rank

    integer :: mpi_err
    
    ! Gather and share the buffers

    call MPI_GATHERV(send_buffer, sendcount, $BUFFER_MPI_TYPE, recv_buffer, recvcounts, displs, &
                     $BUFFER_MPI_TYPE, root_rank, MPI_COMM_WORLD, mpi_err)
    $ASSERT_DEBUG(mpi_err == MPI_SUCCESS,MPI_BCAST failed)

    ! Finish

    return

  end subroutine gatherv_${INFIX}_${BUFFER_RANK}_

  $endsub

  $GATHERV(i_i4,integer(I4),MPI_INTEGER4,0)
  $GATHERV(i_i4,integer(I4),MPI_INTEGER4,1)
  $GATHERV(i_i4,integer(I4),MPI_INTEGER4,2)
  $GATHERV(i_i4,integer(I4),MPI_INTEGER4,3)
  $GATHERV(i_i4,integer(I4),MPI_INTEGER4,4)

  $GATHERV(i_i8,integer(I8),MPI_INTEGER8,0)
  $GATHERV(i_i8,integer(I8),MPI_INTEGER8,1)
  $GATHERV(i_i8,integer(I8),MPI_INTEGER8,2)
  $GATHERV(i_i8,integer(I8),MPI_INTEGER8,3)
  $GATHERV(i_i8,integer(I8),MPI_INTEGER8,4)

  $GATHERV(r_sp,real(SP),MPI_REAL,0)
  $GATHERV(r_sp,real(SP),MPI_REAL,1)
  $GATHERV(r_sp,real(SP),MPI_REAL,2)
  $GATHERV(r_sp,real(SP),MPI_REAL,3)
  $GATHERV(r_sp,real(SP),MPI_REAL,4)

  $GATHERV(r_dp,real(DP),MPI_DOUBLE_PRECISION,0)
  $GATHERV(r_dp,real(DP),MPI_DOUBLE_PRECISION,1)
  $GATHERV(r_dp,real(DP),MPI_DOUBLE_PRECISION,2)
  $GATHERV(r_dp,real(DP),MPI_DOUBLE_PRECISION,3)
  $GATHERV(r_dp,real(DP),MPI_DOUBLE_PRECISION,4)

  $GATHERV(c_sp,complex(SP),MPI_COMPLEX,0)
  $GATHERV(c_sp,complex(SP),MPI_COMPLEX,1)
  $GATHERV(c_sp,complex(SP),MPI_COMPLEX,2)
  $GATHERV(c_sp,complex(SP),MPI_COMPLEX,3)
  $GATHERV(c_sp,complex(SP),MPI_COMPLEX,4)

  $GATHERV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,0)
  $GATHERV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,1)
  $GATHERV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,2)
  $GATHERV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,3)
  $GATHERV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,4)

!****

  $define $ALLGATHERV $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_MPI_TYPE $3
  $local $BUFFER_RANK $4

  subroutine allgatherv_${INFIX}_${BUFFER_RANK}_ (buffer, recvcounts, displs)

    $BUFFER_TYPE, intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: recvcounts(:)
    integer, intent(in)         :: displs(:)

    integer :: mpi_err
    
    ! Gather and share the buffers

    call MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, buffer, recvcounts, displs, &
         $BUFFER_MPI_TYPE, MPI_COMM_WORLD, mpi_err)
    $ASSERT_DEBUG(mpi_err == MPI_SUCCESS,MPI_BCAST failed)

    ! Finish

    return

  end subroutine allgatherv_${INFIX}_${BUFFER_RANK}_

  $endsub

  $ALLGATHERV(i_i4,integer(I4),MPI_INTEGER4,0)
  $ALLGATHERV(i_i4,integer(I4),MPI_INTEGER4,1)
  $ALLGATHERV(i_i4,integer(I4),MPI_INTEGER4,2)
  $ALLGATHERV(i_i4,integer(I4),MPI_INTEGER4,3)
  $ALLGATHERV(i_i4,integer(I4),MPI_INTEGER4,4)

  $ALLGATHERV(i_i8,integer(I8),MPI_INTEGER8,0)
  $ALLGATHERV(i_i8,integer(I8),MPI_INTEGER8,1)
  $ALLGATHERV(i_i8,integer(I8),MPI_INTEGER8,2)
  $ALLGATHERV(i_i8,integer(I8),MPI_INTEGER8,3)
  $ALLGATHERV(i_i8,integer(I8),MPI_INTEGER8,4)

  $ALLGATHERV(r_sp,real(SP),MPI_REAL,0)
  $ALLGATHERV(r_sp,real(SP),MPI_REAL,1)
  $ALLGATHERV(r_sp,real(SP),MPI_REAL,2)
  $ALLGATHERV(r_sp,real(SP),MPI_REAL,3)
  $ALLGATHERV(r_sp,real(SP),MPI_REAL,4)

  $ALLGATHERV(r_dp,real(DP),MPI_DOUBLE_PRECISION,0)
  $ALLGATHERV(r_dp,real(DP),MPI_DOUBLE_PRECISION,1)
  $ALLGATHERV(r_dp,real(DP),MPI_DOUBLE_PRECISION,2)
  $ALLGATHERV(r_dp,real(DP),MPI_DOUBLE_PRECISION,3)
  $ALLGATHERV(r_dp,real(DP),MPI_DOUBLE_PRECISION,4)

  $ALLGATHERV(c_sp,complex(SP),MPI_COMPLEX,0)
  $ALLGATHERV(c_sp,complex(SP),MPI_COMPLEX,1)
  $ALLGATHERV(c_sp,complex(SP),MPI_COMPLEX,2)
  $ALLGATHERV(c_sp,complex(SP),MPI_COMPLEX,3)
  $ALLGATHERV(c_sp,complex(SP),MPI_COMPLEX,4)

  $ALLGATHERV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,0)
  $ALLGATHERV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,1)
  $ALLGATHERV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,2)
  $ALLGATHERV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,3)
  $ALLGATHERV(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,4)

!****

  $define $ALLREDUCE $sub

  $local $INFIX $1
  $local $BUFFER_TYPE $2
  $local $BUFFER_MPI_TYPE $3
  $local $BUFFER_RANK $4

  subroutine allreduce_${INFIX}_${BUFFER_RANK}_ (buffer, op)

    $BUFFER_TYPE, intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
    integer, intent(in)         :: op

    integer :: mpi_err

    ! Reduce the buffers

    call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, PRODUCT(SHAPE(buffer)), $BUFFER_MPI_TYPE, &
         op, MPI_COMM_WORLD, mpi_err)
    $ASSERT_DEBUG(mpi_err == MPI_SUCCESS,MPI_BCAST failed)

    ! Finish

    return

  end subroutine allreduce_${INFIX}_${BUFFER_RANK}_

  $endsub

  $ALLREDUCE(i_i4,integer(I4),MPI_INTEGER4,0)
  $ALLREDUCE(i_i4,integer(I4),MPI_INTEGER4,1)
  $ALLREDUCE(i_i4,integer(I4),MPI_INTEGER4,2)
  $ALLREDUCE(i_i4,integer(I4),MPI_INTEGER4,3)
  $ALLREDUCE(i_i4,integer(I4),MPI_INTEGER4,4)

  $ALLREDUCE(i_i8,integer(I8),MPI_INTEGER8,0)
  $ALLREDUCE(i_i8,integer(I8),MPI_INTEGER8,1)
  $ALLREDUCE(i_i8,integer(I8),MPI_INTEGER8,2)
  $ALLREDUCE(i_i8,integer(I8),MPI_INTEGER8,3)
  $ALLREDUCE(i_i8,integer(I8),MPI_INTEGER8,4)

  $ALLREDUCE(r_sp,real(SP),MPI_REAL,0)
  $ALLREDUCE(r_sp,real(SP),MPI_REAL,1)
  $ALLREDUCE(r_sp,real(SP),MPI_REAL,2)
  $ALLREDUCE(r_sp,real(SP),MPI_REAL,3)
  $ALLREDUCE(r_sp,real(SP),MPI_REAL,4)

  $ALLREDUCE(r_dp,real(DP),MPI_DOUBLE_PRECISION,0)
  $ALLREDUCE(r_dp,real(DP),MPI_DOUBLE_PRECISION,1)
  $ALLREDUCE(r_dp,real(DP),MPI_DOUBLE_PRECISION,2)
  $ALLREDUCE(r_dp,real(DP),MPI_DOUBLE_PRECISION,3)
  $ALLREDUCE(r_dp,real(DP),MPI_DOUBLE_PRECISION,4)

  $ALLREDUCE(c_sp,complex(SP),MPI_COMPLEX,0)
  $ALLREDUCE(c_sp,complex(SP),MPI_COMPLEX,1)
  $ALLREDUCE(c_sp,complex(SP),MPI_COMPLEX,2)
  $ALLREDUCE(c_sp,complex(SP),MPI_COMPLEX,3)
  $ALLREDUCE(c_sp,complex(SP),MPI_COMPLEX,4)

  $ALLREDUCE(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,0)
  $ALLREDUCE(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,1)
  $ALLREDUCE(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,2)
  $ALLREDUCE(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,3)
  $ALLREDUCE(c_dp,complex(DP),MPI_DOUBLE_COMPLEX,4)

  $endif

!****

  subroutine partition_tasks (n, m, k_part)

    integer, intent(in)  :: n
    integer, intent(in)  :: m
    integer, intent(out) :: k_part(:)

    integer :: n_work
    integer :: n_part(SIZE(k_part)-1)
    integer :: i

    n_work = SIZE(k_part) - 1

    ! Partition n tasks among n_work workers --- where possible, at
    ! least m tasks per worker. The resulting partitioning indices are
    ! returned in the array k_part

    n_part = m*(n/(m*n_work))

    size_loop : do i = 1, n_work
       if(SUM(n_part) >= n-m+1) then
          n_part(i) = n_part(i) + n - SUM(n_part)
          exit size_loop
       else
          n_part(i) = n_part(i) + m
       endif
    end do size_loop

    $ASSERT(SUM(n_part) == n,Partitioning failed)

    k_part(1) = 1
      
    index_loop : do i = 1, n_work
       k_part(i+1) = k_part(i) + n_part(i)
    end do index_loop
    
    ! Finish

    return

  end subroutine partition_tasks

end module core_parallel
