.. _troubleshoot-miss:

Missing Modes
=============

For adiabatic oscillation calculations using :program:`gyre`, the
radial order :math:`\numpg` of modes found should be
monotonic-increasing\ [#dipole]_. Departures from this behavior can
occur for a number of reasons.

Insufficient Frequency Resolution
---------------------------------

.. nml:group:: scan
   :no-target:

If the :ref:`frequency grid <freq-grids>` has insufficient resolution,
then :program:`gyre` can skip modes during the bracketing phase, as
discussed in the :ref:`numerical-limits` section. The signature of
insufficient frequency resolution is that an even number of
consecutive modes is missed --- most often, an adjacent pair of modes.

To fix this problem, first check that the distribution of points in
the frequency grids matches (approximately) the expected distribution of
mode eigenfrequencies:

* In the asymptotic limit of large radial order, p modes are uniformly
  distributed in frequency (see, e.g.,
  :ads_citealp:`aerts:2010`). Hence, to search for these modes set
  :nml:option:`grid_type` = :nml:value:`'LINEAR'` in
  the :nml:group:`scan` namelist group(s).

* Likewise, in the asymptotic limit of large radial order, g modes are
  uniformly distributed in period. Hence, to search for these modes
  set :nml:option:`grid_type` = :nml:value:`'INVERSE'`.

* For rotating stars, the asymptotic behaviors mentioned apply in the
  co-rotating reference frame, not in the inertial reference
  frame. So, be sure to also set :nml:option:`grid_frame` =
  :nml:value:`'COROT_I'` or :nml:value:`'COROT_O'`.

Next, try increasing the number of points in the frequency grids,
simply by increasing the :nml:option:`n_freq`.

.. tip::

   A good rule of thumb is that :nml:option:`n_freq` should be around
   5 times larger than the number of modes expected to be found.

Insufficient Spatial Resolution
-------------------------------

.. nml:group:: grid
   :no-target:

If the :ref:`spatial grid <spatial-grids>` has insufficient resolution,
then certain modes can simply be absent from the (finite) set of
distinct numerical solutions, as discussed in the
:ref:`numerical-limits` section. The signature of insufficient spatial
resolution is that modes that `are` found have radial orders
comparable to the number of grid points :math:`N` in the grid; and
that the eigenfunctions of these modes are barely resolved
(cf. :numref:`fig-eigenfuncs-N7`).

To fix this problem, first check that the :nml:option:`w_osc`,
:nml:option:`w_exp`, and :nml:option:`w_ctr` weights are set to
reasonable values (see the :ref:`spatial-grids-rec` section). If that
doesn't improve things, try gradually increasing both :nml:option:`w_osc`
and :nml:option:`w_ctr`.

Non-adiabatic Effects
---------------------

When undertaking :ref:`non-adiabatic calculations <non-ad-osc>`,
modes can be mis-classified or completely missed. The former situation
arises because the expectation of monotonic-increasing :math:`\numpg`
formally applies only to adiabatic oscillations; while it can also
work reasonably well for weakly non-adiabatic cases, there are no
guarantees. If mis-classification does occur, then it must be fixed
manually by determining which adiabatic mode the problematic
non-adiabatic mode corresponds to.

Missing modes occur for a different reason: if a mode has a large
growth rate, then the usual :ref:`adiabatic method <non-ad-ad>`
for establishing initial trial roots can fail to find it. In such
cases, the alternative :ref:`contour method <non-ad-contour>` performs
very well.

.. rubric:: Footnotes

.. [#dipole] The sole exception is :math:`\ell=1` modes, where
             :math:`\numpg=0` is skipped due to the way the
             :ads_citealt:`takata:2006b` classification scheme is set
             up.
