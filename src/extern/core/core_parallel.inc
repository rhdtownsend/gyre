! Incfile  : core_parallel
! Purpose  : parallel support fpx3 macros

$include 'core.inc'

$define $BCAST_1 $sub

$local $BUFFER_TYPE $1

subroutine bcast_1_ (buffer, root_rank)

  $BUFFER_TYPE, intent(inout) :: buffer(:)
  integer, intent(in)         :: root_rank

  integer :: i
  integer :: lb(1)
  integer :: ub(1)

  ! Broadcast the buffer

  lb = LBOUND(buffer)
  ub = UBOUND(buffer)

  do i = lb(1), ub(1)
     call bcast(buffer(i), root_rank)
  end do

  ! Finish

  return

end subroutine bcast_1_

$endsub

!****

$define $BCAST_2 $sub

$local $BUFFER_TYPE $1

subroutine bcast_2_ (buffer, root_rank)

  $BUFFER_TYPE, intent(inout) :: buffer(:,:)
  integer, intent(in)         :: root_rank

  integer :: i
  integer :: j
  integer :: lb(2)
  integer :: ub(2)

  ! Broadcast the buffer

  lb = LBOUND(buffer)
  ub = UBOUND(buffer)

  do j = lb(2), ub(2)
     do i = lb(1), ub(1)
        call bcast(buffer(i,j), root_rank)
     end do
  end do

  ! Finish

  return

end subroutine bcast_2_

$endsub

!****

$define $BCAST_3 $sub

$local $BUFFER_TYPE $1

subroutine bcast_3_ (buffer, root_rank)

  $BUFFER_TYPE, intent(inout) :: buffer(:,:,:)
  integer, intent(in)         :: root_rank

  integer :: i
  integer :: j
  integer :: k
  integer :: lb(3)
  integer :: ub(3)

  ! Broadcast the buffer

  lb = LBOUND(buffer)
  ub = UBOUND(buffer)

  do k = lb(3), ub(3)
     do j = lb(2), ub(2)
        do i = lb(1), ub(1)
           call bcast(buffer(i,j,k), root_rank)
        end do
     end do
  end do

  ! Finish

  return

end subroutine bcast_3_

$endsub

!****

$define $BCAST_4 $sub

$local $BUFFER_TYPE $1

subroutine bcast_4_ (buffer, root_rank)

  $BUFFER_TYPE, intent(inout) :: buffer(:,:,:,:)
  integer, intent(in)         :: root_rank

  integer :: i
  integer :: j
  integer :: k
  integer :: l
  integer :: lb(4)
  integer :: ub(4)

  ! Broadcast the buffer

  lb = LBOUND(buffer)
  ub = UBOUND(buffer)

  do l = lb(4), ub(4)
     do k = lb(3), ub(3)
        do j = lb(2), ub(2)
           do i = lb(1), ub(1)
              call bcast(buffer(i,j,k,l), root_rank)
           end do
        end do
     end do
  end do

  ! Finish

  return

end subroutine bcast_4_

$endsub

!****

$define $BCAST $sub

$local $BUFFER_TYPE $1
$local $BUFFER_RANK $2

$if ($BUFFER_RANK eq 1)
$BCAST_1($BUFFER_TYPE)
$elsif ($BUFFER_RANK eq 2)
$BCAST_2($BUFFER_TYPE)
$elsif ($BUFFER_RANK eq 3)
$BCAST_3($BUFFER_TYPE)
$elsif ($BUFFER_RANK eq 4)
$BCAST_4($BUFFER_TYPE)
$else
$error(Invalid BUFFER_RANK in BCAST)
$endif

$endsub

!****

$define $BCAST_ALLOC $sub

$local $BUFFER_TYPE $1
$local $BUFFER_RANK $2

subroutine bcast_alloc_${BUFFER_RANK}_ (buffer, root_rank)

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

end subroutine bcast_alloc_${BUFFER_RANK}_

$endsub
