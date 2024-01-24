#! Include : forum
#! Purpose : common fypp macros used throughout GYRE

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

#!!
#! Abort execution
#!!

#:def ABORT(message=None)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line ${_LINE_}$ of ${_FILE_}$'
   #:if message
      write(UNIT=ERROR_UNIT, FMT=*) ${message}$
   #:endif
   error stop
#:enddef

#!!
#! Assert that cond is satisfied
#!!

#:def ASSERT(cond, message=None)
   if (.NOT. (${cond}$)) then
      #:set cond_str = cond.replace("'", "''")
      #:if message
         @:ABORT('assertion ${cond_str}$ failed with message '//${message}$//'"')
      #:else
         @:ABORT('assertion ${cond_str}$ failed')
      #:endif
   end if
#:enddef

#!!
#! Same as ASSERT, but only if DEBUG is set
#!!

#:def ASSERT_DEBUG(cond, message=None)
   #:if defined('DEBUG')
      @:ASSERT(${cond}$, ${message}$)
   #:endif
#:enddef

#!!
#! Used for bounds checking
#!!

#:def CHECK_BOUNDS(x, y)
   @:ASSERT_DEBUG(${x}$ == ${y}$, 'bounds mismatch')
#:enddef

#!!
#! Create a deferred shape array specification list
#!!

#:def ARRAY_SPEC(rank, str=':')
   #:if rank > 0
      $:'(' + ','.join([str.format(i) for i in range(1, rank+1)]) + ')'
   #:endif
#:enddef
