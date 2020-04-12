.. _spatial-grids:

Spatial Grids
=============

GYRE discretizes the oscillation equations on a spatial grid
:math:`\{x_{1},x_{2},\ldots,x_{N}\}` in the dimensionless radial
coordinate :math:`x \equiv r/R`. The computational cost of a
calculation scales with the total number of points :math:`N` in this
grid, while the grid's resolution --- i.e., the spacing between
adjacent points --- impacts both the number of modes that can be found
by GYRE, and the accuracy of these modes (see the
:ref:`numerical-limits` section for a discussion of these behaviors in
the context of the stretched string BVP).

Scaffold Grid
-------------

GYRE constructs a fresh spatial grid for each combination of harmonic
degree :math:`\ell` and azimuthal order :math:`m` specified in the
:nml_g:`mode` namelist groups (see the :ref:`input-files` chapter for
more details). The starting point for each of these grids is the
*scaffold grid*, which comprises the following:

* an inner point :math:`x=\xin`;
* an outer point :math:`x=\xout`;
* the subset of points of the input model grid satisfying :math:`\xin <
  x < \xout`

By default, :math:`\xin` and :math:`\xout` are obtained from the input
model grid as well, meaning that the scaffold grid is identical to the
model grid. However, either or both can be overridden using the
:nml_n:`x_i` and :nml_n:`x_o` parameters, respectively, of the
:nml_g:`grid` namelist group.

Iterative Refinement
--------------------

GYRE refines a scaffold grid through a sequence of iterations. During
a given iteration, each subinterval :math:`[x_{k},x_{k+1}]`
(:math:`k=1,2,\ldots,N-1`) is assessed against various criteria
(discusssed in greater detail below). If any criteria match, then the
subinterval is refined by bisection, inserting an additional point at
the midpoint

.. math::

   x_{k+1/2} = \frac{x_{k} + x_{k+1}}{2}.

The sequence terminates if no refinements occur during a given
iteration, or if the number of completed iterations equals the value
specified by the :nml_n:`n_iter_max` parameter of the :nml_g:`grid`
namelist group.

.. _wave-criterion:

Mechanical Criterion
~~~~~~~~~~~~~~~~~~~~

The wave criterion involves a local analysis of the mechanical parts
of the oscillation equations, with the goal of improving resolution
where the displacement perturbation :math:`\vxi` is rapidly
varying. Within the subinterval :math:`[x_{k},x_{k+1}]`, the
:math:`y_{1}` and :math:`y_{2}` solutions (see the
:ref:`math-formalism` chapter) take the approximate form

.. math::

   y_{1,2}(x) \sim \exp [ \chi \, (\ln x - \ln x_{k+1/2}) ],

where :math:`\chi` is one of the eigenvalues of the mechanical
(upper-left) :math:`2 \times 2` submatrix of the full Jacobian matrix
:math:`\mA` , evaluated at the midpoint :math:`x_{k+1/2}`.

In propagation zones the imaginary part :math:`\chi_{\rm i}` of the
eigenvalue gives the local wavenumber in :math:`\ln x` space, and
:math:`2\pi \chi_{\rm i}^{-1}` the corresponding wavelength; while in
evanescent zones the real part :math:`\chi_{\rm r}` gives the local
exponential growth/decay rate, and :math:`\chi_{\rm r}^{-1}` the
corresponding e-folding length.

Based on this analysis, the criterion for refinement of the
subinterval is

.. math::

   ( \ln x_{k+1} - \ln x_{k} ) \, \max (\alpha_{\rm osc} |\chi_{\rm i}|, \alpha_{\rm exp} |\chi_{\rm r}|) > 2 \pi

This causes refinement if the subinterval width (in :math:`\ln x`
space) exceeds :math:`\alpha_{\rm osc}^{-1}` times the local
wavelength, or :math:`2\pi \alpha_{\rm exp}^{-1}` times the local
e-folding length. The controls :math:`\alpha_{\rm exp}` and
:math:`\alpha_{\rm exp}` are set via the :nml_n:`alpha_exp` and
:nml_n:`alpha_osc` parameters, respectively, of the :nml_g:`grid`
namelist group.

.. tip::

   While :nml_n:`alpha_exp` and :nml_n:`alpha_osc` default to zero, it
   is highly recommended to use non-zero values for these parameters,
   to ensure adequate resolution of solutions throughout the
   star. Reasonable starting choices are :nml_n:`alpha_osc = 10` and
   :nml_nv:`alpha_exp = 2`.

Because there are two possible values for :math:`\chi`, the above
refinement criterion is applied twice (once for each). Moreover,
because :math:`\chi` depends implicitly on the oscillation frequency,
the criterion is applied for each frequency in the grid
:math:`\{\omega_{1},\omega_{2},\ldots,\omega_{M}\}`.

.. _thermal-criterion:

Thermal Criterion
~~~~~~~~~~~~~~~~~

Similar to the wave criterion discussed above, the thermal criterion
involves a local analysis of the energetic parts of the oscillation
equation, with the goal of improving resolution where the thermal
timescale is very long and perturbations are almost adiabatic. Within
the subinterval :math:`[x_{k},x_{k+1}]`, the :math:`y_{5}` and
:math:`y_{6}` perturbation take the approximate form

.. math::

   y_{5,6}(x) \sim \exp [ \pm \tau \, (\ln x - \ln x_{k+1/2}) ],

where :math:`\pm\tau` are the eigenvalues of the matrix formed from
the energetic (bottom-rright) :math:`2 \times 2` submatrix of the full
Jacobian matrix :math:`\mA`, evaluated at the midpoint
:math:`x_{k+1/2}`.

Based on this analysis, the criterion for refinement of the
subinterval is

.. math::

   ( \ln x_{k+1} - \ln x_{k} ) \, \alpha_{\rm thm} |\tau| > 1.

The control :math:`\alpha_{\rm thm}` is set via the :nml_n:`alpha_thm`
parameters of the :nml_g:`grid` namelist group.

Because :math:`\tau` depends implicitly on the oscillation frequency,
this criterion is applied for each frequency in the grid
:math:`\{\omega_{1},\omega_{2},\ldots,\omega_{M}\}`.

.. _structural-criteria:

Structural Criteria
~~~~~~~~~~~~~~~~~~~

The structural criteria have the goal of improving resolution where
the stellar structure coefficients are changing rapidly. For a given
coefficient :math:`C`, the criterion for refinement of the subinterval
:math:`[x_{k},x_{k+1}]` is

.. math::

   ( \ln x_{k+1} - \ln x_{k} ) \, \alpha_{\rm str} \left| \pderiv{\ln C}{\ln x} \right| > 1

The control :math:`\alpha_{\rm thm}` is set via the :nml_n:`alpha_thm`
parameter of the :nml_g:`grid` namelist group. This criterion is
applied to the :math:`V_2 \equiv V/x`, :math:`U`, :math:`A^{*}`,
:math:`c_{1}` and :math:`\Gamma_{1}` coefficients (see the
:ref:`structure-coeffs` section).

.. _central-criteria:

Central Criteria
~~~~~~~~~~~~~~~~

All of the above criteria depend on the logarithmic subinterval width
:math:`(\ln x_{k+1} - \ln x_{k})`, and cannot be applied to the first
subinterval :math:`[x_{1},x_{2}]` if it extends to the center of the
star :math:`x = 0`. In such cases, the :nml_n:`resolve_ctr` parameter
of the :nml_g:`grid` namelist group determines whether the subinterval
is refined. If set to :nml_v:`.FALSE.`, then no refinement occurs;
while if set to :nml_v:`.TRUE.`, then the refinement criteria are

.. math::

   \chi_{\rm i} > 0

or

.. math::

   \alpha_{\rm ctr} | \chi_{rm r} | > 1

where :math:`\chi` is the eigenvalue from the local analysis (see the
:ref:`wave-criterion` section) corresponding to the solution that
remains well-behaved at the origin. The first criterion causes
refinement if the subinterval is in a propagation zone, and the second
if the solution slope :math:`|\sderiv{\ln y}{\ln x}| \sim |\chi_{\rm
r}|` exceeds :math:`\alpha_{\rm ctr}^{-1}`. The control
:math:`\alpha_{\rm ctr}` is set via the :nml_n:`alpha_ctr` parameters
of the :nml_g:`grid` namelist group.

.. tip::

   While :nml_n:`alpha_ctr` defaults to zero, it is highly recommended
   to use a non-zero value for this parameter, to ensure adequate
   resolution of solutions at the center. A reasonable starting choice
   is :nml_n:`alpha_ctr = 10`.

Because :math:`\chi` depends implicitly on the oscillation frequency,
these criteria are applied for each frequency in the grid
:math:`\{\omega_{1},\omega_{2},\ldots,\omega_{M}\}`.
