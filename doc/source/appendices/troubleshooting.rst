.. _troubleshooting:

***************
Troubleshooting
***************

Missng Modes
============

For adiabatic calculations the radial order :math:`\npg` of modes
found, as GYRE processes a given :nml_g:`mode` namelist group, should
be monotonic-increasing\ [#dipole]_. Departures from this behavior can
occur for a number of reasons:

  * Missing values can indicate that GYRE has skipped a mode in
    frequency space; the fix is to use a finer frequency grid.

  * Missing values together with duplicate and/or non-monotonic values
    can indicate that GYRE isn't resolving the spatial structure of
    eigenfunctions; the fix is to use a finer spatial grid.

  * Missing values together with duplicate and/or non-monotonic values
    can *also* incdicate problems with the input stellar model ---
    for instance, incorrect values for the Brunt-Vaisala frequency
    across density discontinuities; the fix is to stop expecting GYRE
    to give sensible output when fed crap stellar models!

.. rubric:: Footnotes

.. [#dipole] The sole exception is :math:`\ell=1` modes, where
             :math:`\npg=0` is skipped due to the way the
             :ads_citealt:`takata:2006b` classification scheme is set
             up.
