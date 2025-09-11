#! Include : gyre
#! Purpose : common macros used throughout GYRE
#!
#! Copyright 2013-2025 Rich Townsend & The GYRE Team
#!
#! This file is part of GYRE. GYRE is free software: you can
#! redistribute it and/or modify it under the terms of the GNU General
#! Public License as published by the Free Software Foundation, version 3.
#!
#! GYRE is distributed in the hope that it will be useful, but WITHOUT
#! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
#! License for more details.
#!
#! You should have received a copy of the GNU General Public License
#! along with this program.  If not, see <http://www.gnu.org/licenses/>.

#!!! Include other files

#:include 'forum.inc'

#!!! Type and suffix definitions

#:set INT_TYPES = ['integer']
#:set INT_SUFFIXES = ['i']

#:set REAL_TYPES = ['real(RD)']
#:set REAL_SUFFIXES = ['r']

#:set COMPLEX_TYPES = ['complex(RD)']
#:set COMPLEX_SUFFIXES = ['c']

#:set FLOAT_TYPES = REAL_TYPES + COMPLEX_TYPES
#:set FLOAT_SUFFIXES = REAL_SUFFIXES + COMPLEX_SUFFIXES

#:set EXT_REAL_TYPES = ['type(ext_rt)']
#:set EXT_REAL_SUFFIXES = ['xr']

#:set EXT_COMPLEX_TYPES = ['type(ext_ct)']
#:set EXT_COMPLEX_SUFFIXES = ['xc']

#:set EXT_FLOAT_TYPES = EXT_REAL_TYPES + EXT_COMPLEX_TYPES
#:set EXT_FLOAT_SUFFIXES = EXT_REAL_SUFFIXES + EXT_COMPLEX_SUFFIXES

#:set CHAR_TYPES = ['character(*)']
#:set CHAR_SUFFIXES = ['a']

#:set LOGICAL_TYPES = ['logical']
#:set LOGICAL_SUFFIXES = ['l']

#:set NUM_TYPES = INT_TYPES + REAL_TYPES + COMPLEX_TYPES
#:set NUM_SUFFIXES = INT_SUFFIXES + REAL_SUFFIXES + COMPLEX_SUFFIXES

#:set MATH_TYPES = NUM_TYPES + LOGICAL_TYPES
#:set MATH_SUFFIXES = NUM_SUFFIXES + LOGICAL_SUFFIXES

#:set DTYPE_BASE_TYPES = ['real(RD)', 'complex(RD)']
#:set DTYPE_SUFFIXES = ['rt', 'ct']

#:set PAR_NAMES = ['grid', 'num', 'orbit', 'osc', 'rot', 'scan']
#:set PAR_TYPES = [f'type({name}_par_t)' for name in PAR_NAMES]
#:set PAR_SUFFIXES = ['gr', 'nm', 'or', 'os', 'rt', 'sc']

#!!! Macro definitions

#!!
#! Print out a warning
#!!

#:def WARN(message=None)
  #:if message
    write(UNIT=ERROR_UNIT, FMT=*) 'INFO at line ${_LINE_}$ of ${_FILE_}$'
    write(UNIT=ERROR_UNIT, FMT=*) ${message}$
  #:else
    write(UNIT=ERROR_UNIT, FMT=*) 'INFO at line ${_LINE_}$ of ${_FILE_}$'
  #:endif
#:enddef
