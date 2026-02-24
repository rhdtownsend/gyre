.. _num-group:

.. nml:group:: num

Num Namelist Group
==================

The :nml:group:`num` namelist group controls numerical aspects of
calculations. The input file can contain one or more, but only the
last (tag-matching) one is used. The following options are available:

.. nml:option:: diff_scheme
   :type: string
   :default: 'COLLOC_GL2'`

   Difference equation scheme; one of

   - :nml:value:`'COLLOC_GL2'` : Second-order Gauss-Legendre collocation
   - :nml:value:`'COLLOC_GL4'` : Fourth-order Gauss-Legendre collocation
   - :nml:value:`'COLLOC_GL6'` : Sixth-order Gauss-Legendre collocation
   - :nml:value:`'MAGNUS_GL2'` : Second-order Gauss-Legendre Magnus
   - :nml:value:`'MAGNUS_GL4'` : Fourth-order Gauss-Legendre Magnus
   - :nml:value:`'MAGNUS_GL6'` : Sixth-order Gauss-Legendre Magnus
   - :nml:value:`'MIRK'` : Fourth-order mono-implicit Runge-Kutta (experimental)
   - :nml:value:`'TRAPZ'` : Trapezoidal, with the prescription by :ads_citet:`sugimoto:1970` for non-adiabatic cases

.. nml:option:: real_root_solver
   :type: string
   :default: 'BRENT'

   Root solver for real arithmetic; one of

   - :nml:value:`'BRENT'` : Brent's method

.. nml:option:: cmplx_root_solver
   :type: string
   :default: 'RIDDERS'

   Root solver for complex arithmetic; one of

   - :nml:value:`'RIDDERS'` : Complex Ridders' method
   - :nml:value:`'SECANT'` : Secant method
   - :nml:value:`'SIMPLEX'` : Simplex method

.. nml:option:: n_iter_max
   :type: integer
   :default: 50

   Maximum number of iterations in root-finding algorithms

.. nml:option:: ad_matrix_solver
   :type: string
   :default: 'ROWPP'`

   Matrix solver for discretized adiabatic equations; one of

   - :nml:value:`'BANDED'` : Banded factorization (LAPACK xGBTRF routines)
   - :nml:value:`'CYCLIC'` : Cyclic structured factorization (`Wright 1994 <https://link.springer.com/article/10.1007/s002110050043>`__)
   - :nml:value:`'ROWPP'`  : Gaussian elimination with row partial pivoting

.. nml:option:: nad_matrix_solver
   :type: string
   :default: 'CYCLIC'

   Matrix solver for discretized non-adiabatic equations; same choices
   as for :nml:option:`ad_matrix_solver`

.. nml:option:: deflate_roots
   :type: logical
   :default: .TRUE.

   Use deflation during complex root finding; this helps avoid the
   same eigenfrequency being found multiple times

.. nml:option:: restrict_roots
   :type: logical
   :default: .TRUE.

   Discard roots that fall outside the bounds of the frequency scan

.. nml:option:: ad_search
   :type: string
   :default: 'BRACKET'

   Initial search method for adiabatic calculations; one of

   - :nml:value:`'BRACKET'` : Bracket sign changes in the discriminant function

.. nml:option:: nad_search
   :type: string
   :default: 'AD'

   Initial search method for non-adiabatic calculations; one of

   - :nml:value:`'AD'` : Use adiabatic eigenfrequencies
   - :nml:value:`'MINMOD'` : Find minima in the modulus of the discriminant function, along the real-:math:`\omega` axis
   - :nml:value:`'CONTOUR'` : Find intersections between real and imaginary zero-contours of the discriminant function

   See the :ref:`non-ad-osc` chapter for further details about these search methods

.. nml:option:: tag_list
   :type: string
   :default: ''

   Comma-separated list of :nml:option:`tag <mode.tag>` values to match;
   matches all if left blank
