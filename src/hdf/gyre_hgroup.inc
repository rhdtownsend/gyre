! Incfile  : gyre_hgroup
! Purpose  : macros for gyre_hgroup module
!
! Copyright 2020 Rich Townsend
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

$perl
xmacro 'HDF5_CALL', sub {
  my @a = @_; @a = get_args() unless @a;
  die("Invalid number of arguments:".scalar @a) unless scalar @a >= 1;
  my $proc_name = shift @a;
  @pos_args = grep(!/=/, @a);
  @key_args = grep(/=/, @a);
  @args = (@pos_args, 'hdf_err', @key_args);
  my $arg_list = join(', ', @args);
  return <<EOF;
call $proc_name($arg_list)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line $line_num $context_doc:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to $proc_name'
   stop 'Program aborted'
endif
EOF
};
$endperl
