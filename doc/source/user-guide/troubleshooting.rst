.. _troubleshooting:

***************
Troubleshooting
***************

This chapter discusses various problems that can arise during normal
GYRE operation, and steps that can be taken to resolve them.

Missing Modes
=============

For adiabatic calculations the radial order :math:`\npg` of modes
found, as GYRE processes a given :nml_g:`mode` namelist group, should
be monotonic-increasing\ [#dipole]_. Departures from this behavior can
occur for a number of reasons, as discussed in the following subsections.

Insufficient Frequency Resolution
---------------------------------

If the :ref:`frequency grid <freq-grids>` has insufficient resolution,
then GYRE can skip modes during the bracketing phase, as discussed in
the :ref:`numerical-limits` section. The signature of insufficient
frequency resolution is that an even number of consecutive modes is missed ---
most often, an adjacent pair of modes.

To fix this problem, first check that the distribution of points in
the frequency grids matches (appromately) the expected distribution of
mode eigenfrequencies:

* In the asymptotic limit of large radial order and/or harmonic
  degree, p modes are uniformly distributed in frequency (see, e.g.,
  :ads_citealp:`aerts:2010`). Hence, to search for these modes set
  :nml_n:`grid_type`\ =\ :nml_v:`'LINEAR'` in the :nml_g:`scan`
  namelist group(s).

* Likewise, in the asymptotic limit of large radial order and/or
  harmonic degree, g modes are uniformly distributed in period. Hence,
  to search for these modes set :nml_n:`grid_type`\ =\
  :nml_v:`'INVERSE'` in the :nml_g:`scan` namelist group(s).

* For rotating stars, the asymptotic behaviors mentiond apply in the
  co-rotating reference frame, not in the inertial reference
  frame. So, be sure to also set :nml_n:`grid_frame` \ =\
  :nml_n:`'COROT_I'`\ \|\ :nml_n:`'COROT_O'` in the :nml_g:`scan`
  namelist group.
  

Next, try increasing the number of points in the frequency grids,
simply by increasing the :nml_n:`n_freq` parameter in the
:nml_g:`scan` namelist group(s).

.. tip::

   A good rule of thumb is that :nml_n:`n_freq` should be around 5
   times larger than the number of modes expected to be found.

Insufficient Spatial Resolution
-------------------------------

If the :ref:`spatial grid <freq-grids>` has insufficient resolution,
then certain modes can simply be absent from the (finite) set of
distinct numerical solutions, as discussed in the
:ref:`numerical-limits` section. The signature of insufficient spatial
resolution is that modes that `are` found have radial orders
comparable to the number of grid points :math:`N` in the grid; and
that the eigenfunctions of these modes are barely resolved
(cf. :numref:`fig-eigenfuncs-N7`).

To fix this issue

Duplicated Modes
================



.. rubric:: Footnotes

.. [#dipole] The sole exception is :math:`\ell=1` modes, where
             :math:`\npg=0` is skipped due to the way the
             :ads_citealt:`takata:2006b` classification scheme is set
             up.


	     
