.. _long-runtimes:

Long Runtimes
=============

gyre
----

Long runtimes with the :program:`gyre` frontend occur when the
:ref:`spatial grid <spatial-grids>` and/or :ref:`frequency grid
<freq-grids>` contain many points. The execution time to process a
single :nml_g:`mode` namelist group can be approximated by

.. math::

   \tau \approx C_{\rm b} N M + C_{\rm s} N N_{\rm m},

where :math:`N` is the number of spatial points, :math:`M` is the
number of frequency points, :math:`N_{\rm m}` is the number of modes
found, and :math:`C_{\rm b}` and :math:`C_{\rm s}` are constants. The
first (:math:`C_{\rm b}`) term represents the time take to bracket
roots of the discriminant function, and the second (:math:`C_{\rm s}`)
the time taken to solve for these roots (see the
:ref:`numerical` chapter for details).

The key to ensuring reasonable runtimes lies in judicious choice of
parameters in the :nml_g:`scan` namelist group(s). The :nml_n:`n_freq`
parameter obviously has an impact on :math:`\tau`, as it directly sets
:math:`M`. However, the :nml_n:`freq_min` and :nml_n:`freq_max`
parameters also influence :math:`\tau`, due to the way the spatial
grid is constructed. If the frequency scan includes parts of the
star's oscillation spectrum containing modes with very large radial
orders (whether p modes or g modes), then GYRE's :ref:`iterative
refinement algorithm <spatial-grids-iter>` will insert many grid
points in order to resolve the modes' wavefunctions. This can
ultimately lead to huge :math:`N` and very long runtimes.

gyre_tides
----------

The narrative is similar with the :program:`gyre_tides` frontend,
although there are no frequency grids involved. The execution time to
process a single :nml_g:`orbit` namelist group can be approximated by

.. math::

   \tau \approx C_{\rm t} N

where :math:`C_{\rm t}` is a constant.






