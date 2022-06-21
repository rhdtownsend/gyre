! Incfile  : gyre_detail
! Purpose  : support macros for gyre_detail
!
! Copyright 2021-2022 Rich Townsend & The GYRE Team
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
  
$define $WRITE_POINTS $sub
$local $NAME $1
$local $FUNC $2

case ($str($NAME))

   call wr%write($str($NAME), [($FUNC, j=1,gr%n)])

$endsub

$define $WRITE_VALUE $sub
$local $NAME $1
$local $VALUE $2

case ($str($NAME))

   call wr%write($str($NAME), $VALUE)

$endsub
