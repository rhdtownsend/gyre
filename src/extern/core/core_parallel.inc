! Incfile  : core_parallel
! Purpose  : parallel support fpx3 macros

$include 'core.inc'

$define $BCAST_ALLOC $sub

$local $INFIX $1
$local $BUFFER_TYPE $2
$local $BUFFER_RANK $3

subroutine bcast_alloc_${INFIX}_${BUFFER_RANK} (buffer, root_rank)

  $BUFFER_TYPE, allocatable, intent(inout) :: buffer$ARRAY_SPEC($BUFFER_RANK)
  integer, intent(in)                      :: root_rank

  logical :: alloc
  $if($BUFFER_RANK > 0)
  integer :: lb($BUFFER_RANK)
  integer :: ub($BUFFER_RANK)
  $endif

  ! Deallocate the buffer on non-root processors

  if(MPI_RANK /= root_rank .AND. ALLOCATED(buffer)) then
     deallocate(buffer)
  endif

  ! Check if the buffer is allocated on the root processor

  if(MPI_RANK == root_rank) alloc = ALLOCATED(buffer)
  call bcast(alloc, root_rank)

  if(alloc) then

     ! Broadcast the buffer bounds

     $if($BUFFER_RANK > 0)

     if(MPI_RANK == root_rank) then
        lb = LBOUND(buffer)
        ub = UBOUND(buffer)
     endif

     call bcast(lb, root_rank)
     call bcast(ub, root_rank)

     $endif

     ! Allocate the buffer

     $if($BUFFER_RANK > 0)

     if(MPI_RANK /= root_rank) allocate(buffer($ARRAY_EXPAND(lb,ub,$BUFFER_RANK)))

     $else

     if(MPI_RANK /= root_rank) allocate(buffer)

     $endif

     ! Broadcast the buffer

     call bcast(buffer, root_rank)

  endif

  ! Finish

  return

end subroutine bcast_alloc_${INFIX}_${BUFFER_RANK}

$endsub
