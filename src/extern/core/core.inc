! Incfile  : core
! Purpose  : fpx3 macros

$define $ABORT $sub

  $local $MESSAGE $1

  $if($MESSAGE)
    write(UNIT=ERROR_UNIT, FMT=*) $str(ABORT at line $line_num $context_doc:)
    write(UNIT=ERROR_UNIT, FMT=*) $str($MESSAGE)
  $else
    write(UNIT=ERROR_UNIT, FMT=*) $str(ABORT at line $line_num $context_doc)
  $endif

  stop 'Program aborted'

$endsub

!****

$define $ABORT_DEBUG $sub

  $local $MESSAGE $1

  $if($DEBUG)
    $ABORT($MESSAGE)
  $endif

$endsub

!****

$define $WARN $sub

  $local $MESSAGE $1

  $if($MESSAGE)
    write(UNIT=ERROR_UNIT, FMT=*) $str(INFO at line $line_num $context_doc:)
    write(UNIT=ERROR_UNIT, FMT=*) $str($MESSAGE)
  $else
    write(UNIT=ERROR_UNIT, FMT=*) $str(INFO at line $line_num $context_doc)
  $endif

$endsub

!****

$define $WARN_DEBUG $sub

  $local $MESSAGE $1

  $if($DEBUG)
    $WARN($MESSAGE)
  $endif

$endsub

!****

$define $ASSERT $sub

  $local $CONDITION $1
  $local $MESSAGE $2

  $if($MESSAGE)
    if(.NOT. ($CONDITION)) then
      write(UNIT=ERROR_UNIT, FMT=*) $str(ASSERT ''$CONDITION'' failed at line $line_num $context_doc:)
      write(UNIT=ERROR_UNIT, FMT=*) $str($MESSAGE)
      stop
    endif
  $else
    if(.NOT. $CONDITION) then
      write(UNIT=ERROR_UNIT, FMT=*) $str(ASSERT ''$CONDITION'' failed at line $line_num $context_doc)
      stop
    endif
  $endif

$endsub

!****

$define $ASSERT_DEBUG $sub

  $local $CONDITION $1
  $local $MESSAGE $2

  $if($DEBUG)
    $ASSERT($CONDITION,$MESSAGE)
  $endif

$endsub

!****

$define $CHECK_BOUNDS $sub

  $local $SIZE_A $1
  $local $SIZE_B $2

  $if($DEBUG)

  if($SIZE_A/= $SIZE_B) then
    write(UNIT=ERROR_UNIT, FMT=*) $str($SIZE_A :), $SIZE_A
    write(UNIT=ERROR_UNIT, FMT=*) $str($SIZE_B :), $SIZE_B
    write(UNIT=ERROR_UNIT, FMT=*) $str(CHECK_BOUNDS $SIZE_A==$SIZE_B failed at line $line_num $context_doc)
    stop
  endif

  $endif

$endsub

!****

$perl
xmacro 'ARRAY_SPEC', sub {
  my @a = @_; @a = get_args() unless @a;
  die("Invalid number of arguments:".scalar @a) unless (scalar @a == 1 || scalar @a == 2);
  if($a[0] > 0) {
    my @dim_spec = ();
    for(my $i=1; $i <= $a[0]; $i++) {
      if(defined($a[1])) {
	push @dim_spec, "SIZE($a[1],$i)";
      }
      else {
	push @dim_spec, ":";
      }
    }
    return "(".join(",", @dim_spec).")";
  }
  else {
    return "";
  }
};
$endperl

!****

$perl
xmacro 'ARRAY_EXPAND', sub {
  my @a = @_; @a = get_args() unless @a;
  die("Invalid number of arguments:".scalar @a) unless (scalar @a == 2 || scalar @a == 3);
  my @elems = ();
  if(scalar @a == 3) {
      for(my $i=1; $i <= $a[2]; $i++) {
	  push @elems, "$a[0]($i):$a[1]($i)";
      }
  } else {
      for(my $i=1; $i <= $a[1]; $i++) {
	  push @elems, "$a[0]($i)";
      }
  }
  return join(",", @elems);
};
$endperl

!****

$perl
xmacro 'ARRAY_FIRST', sub {
  my @a = @_; @a = get_args() unless @a;
  die("Invalid number of arguments") if scalar @a != 2;
  if($a[1] > 0) {
      return $a[0]."(".join(",", ("1") x $a[1]).")";
  }
  else {
      return $a[0];
  }
};
$endperl
