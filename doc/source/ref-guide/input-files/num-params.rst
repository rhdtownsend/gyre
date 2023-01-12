.. _num-params:

Numerical Parameters
====================

The :nml_g:`num` namelist group defines numerical method parameters; the
input file can contain one or more, but only the last (tag-matching) one is
used. Allowable fields are:

:nml_n:`diff_scheme` (default :nml_v:`'COLLOC_GL2'`)
  Difference equation scheme; one of:

  - :nml_v:`'COLLOC_GL2'` : Second-order Gauss-Legendre collocation
  - :nml_v:`'COLLOC_GL4'` : Fourth-order Gauss-Legendre collocation
  - :nml_v:`'COLLOC_GL6'` : Sixth-order Gauss-Legendre collocation
  - :nml_v:`'MAGNUS_GL2'` : Second-order Gauss-Legendre Magnus
  - :nml_v:`'MAGNUS_GL4'` : Fourth-order Gauss-Legendre Magnus
  - :nml_v:`'MAGNUS_GL6'` : Sixth-order Gauss-Legendre Magnus
  - :nml_v:`'MIRK'` : Fourth-order mono-implicit Runge-Kutta (experimental)
  - :nml_v:`'TRAPZ'` : Trapezoidal, with the prescription by
    :ads_citet:`sugimoto:1970` for non-adiabatic cases

:nml_n:`r_root_solver` (default :nml_v:`'BRENT'`)
  Root solver for real arithmetic; one of:

  - :nml_v:`'BRENT'` : Brent's method

:nml_n:`c_root_solver` (default :nml_v:`'RIDDERS'`)
  Root solver for complex arithmetic; one of

  - :nml_v:`'RIDDERS'` : Complex Ridders' method
  - :nml_v:`'SECANT'` : Secant method
  - :nml_v:`'SIMPLEX'` : Simplex method

:nml_n:`n_iter_max` (default :nml_v:`50`)
  Maximum number of iterations in root-finding algorithm
  
:nml_n:`matrix_type` (default :nml_v:`'BLOCK`')
  Storage type of system matrix; one of

  - :nml_v:`'BAND'` : Band-structured
  - :nml_v:`'BLOCK'` : Block-structured

:nml_n:`deflate_roots` (default :nml_v:`.TRUE.`)
  Flag to use root deflation, which can avoid the same eigenfrequency
  being found multiple times

:nml_n:`restrict_roots` (default :nml_v:`.TRUE.`)
  Flag to check each roots found lies within the bounds of the frequency scan

:nml_n:`ad_search` (default :nml_v:`'BRACKET'`)
  Initial search method for adiabatic calculations; one of

  - :nml_v:`'BRACKET'` : Bracket sign changes in the discriminant function 

:nml_n:`nad_search` (default :nml_v:`'AD'`)
  Initial search method for non-adiabatic calculations; one of

  - :nml_v:`'AD'` : Use adiabatic eigenfrequencies
  - :nml_v:`'MINMOD'` : Find minima in the modulus of the discriminant function, along the real-:math:`\omega` axis
  - :nml_v:`'CONTOUR'` : Find intersections between real and imaginary zero-contours of the discriminant function

  See the :ref:`non-ad-osc` chapter for more details about these search methods.
    
:nml_n:`tag_list` (default :nml_v:`''`, which matches all)
   Comma-separated list of :nml_g:`mode` tags to match
