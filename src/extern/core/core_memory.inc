! Incfile  : core_memory
! Purpose  : memory management fpx3 macros

$include 'core.inc'

$perl
xmacro 'COPY_ARRAY', sub {
  my @a = @_; @a = get_args() unless @a;
  die("Invalid number of arguments:".scalar @a) unless scalar @a == 1;
  my $rank = shift @a;
  my $output = '';
  for(my $i = 0; $i < 2**$rank; $i++) {
    my @spec_left;
    my @spec_right;
    for(my $j = 0; $j < $rank; $j++) {
      my $k = $j + 1;
      if($i & 2**$j) {
        push @spec_left, "i_c($k):";
        push @spec_right, "i_b($k):";
      }
      else {
        push @spec_left, ":i_a($k)";
        push @spec_right, ":i_a($k)";
      }
    }
    $output .= "array_new(".join(',', @spec_left).") = array(".join(',', @spec_right).")\n";
  }
  return $output;
};
$endperl

$define $REALLOCATE $sub

$if ($sub_arg_num eq 3)
$local $INFIX $1_$3
$local $ARRAY_TYPE $2
$local $ARRAY_RANK $3
$else
$local $INFIX $2
$local $ARRAY_TYPE $1
$local $ARRAY_RANK $2
$endif

$local $ARRAY_NEW_TYPE $eval($ARRAY_TYPE eq 'character(*)' ? 'character(LEN(array))' : $ARRAY_TYPE)

subroutine reallocate_${INFIX}_ (array, shape_new, start)

  $ARRAY_TYPE, allocatable, intent(inout) :: array$ARRAY_SPEC($ARRAY_RANK)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  $ARRAY_NEW_TYPE, allocatable :: array_new$ARRAY_SPEC($ARRAY_RANK)
  integer                      :: i_a($ARRAY_RANK)
  integer                      :: i_b($ARRAY_RANK)
  integer                      :: i_c($ARRAY_RANK)

  $ASSERT_DEBUG(SIZE(shape_new) == $ARRAY_RANK,Dimension mismatch)

  if(PRESENT(start)) then
     $ASSERT_DEBUG(SIZE(start) == $ARRAY_RANK,Dimension mismatch)
  end if
    
  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new($ARRAY_EXPAND(shape_new,$ARRAY_RANK)))

  ! Copy the data from old to new

  $COPY_ARRAY($ARRAY_RANK)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_${INFIX}_

$endsub

