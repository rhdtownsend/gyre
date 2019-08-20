Numerical Parameters
====================

The :nml_g:`num` namelist group defines numerical method parameters; the
input file can contain one or more, but only the last (matching) one is
used. Allowable fields are:

:nml_o:`diff_scheme` (default :nml_l:`'COLLOC_GL2'`)
  Difference equation scheme; one of:

  - :nml_l:`'COLLOC_GL2'` : Second-order Gauss-Legendre collocation
  - :nml_l:`'COLLOC_GL4'` : Fourth-order Gauss-Legendre collocation
  - :nml_l:`'COLLOC_GL6'` : Sixth-order Gauss-Legendre collocation
  - :nml_l:`'MAGNUS_GL2'` : Second-order Gauss-Legendre Magnus
  - :nml_l:`'MAGNUS_GL4'` : Sourth-order Gauss-Legendre Magnus
  - :nml_l:`'MAGNUS_GL6'` : Sixth-order Gauss-Legendre Magnus
  - :nml_l:`'MIRK'` : Fourth-order mono-implicit Runge-Kutta (experimental)
  - :nml_l:`'TRAPZ'` : Trapezoidal, with the :ads:`Sugimoto (1970) <1970ApJ...159..619S>` prescription for non-adiabatic cases

:nml_o:`r_root_solver` (default :nml_l:`'BRENT'`)
  Root solver for real arithmetic; one of:

  - :nml_l:`'BRENT'` : Brent's method

:nml_o:`c_root_solver` (default :nml_l:`'RIDDERS'`)
  Root solver for complex arithmetic; one of

  - :nml_l:`'RIDDERS'` : Complex Ridders' method
  - :nml_l:`'SECANT'` : Secant method
  - :nml_l:`'SIMPLEX'` : Simplex method

:nml_o:`n_iter_max` (default :nml_l:`50`)
  Maximum number of iterations in root-finding algorithm
  
:nml_o:`matrix_type` (default :nml_l:`'BLOCK`')
  Storage type of system matrix; one of

  - :nml_l:`'BAND'` : Band-structured
  - :nml_l:`'BLOCK'` : Block-structued

:nml_o:`deflate_roots` (default :nml_l:`.TRUE.`)
  Flag to use root deflation, which can avoid the same eigenfrequency
  being found multiple times

:nml_o:`restrict_roots` (default :nml_l:`.TRUE.`)
  Flag to check each roots found lies within the bounds of the frequency scan

:nml_o:`tag_list` (default :nml_l:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match
