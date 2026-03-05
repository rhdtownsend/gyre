.. _perfomance:

.. nml:group:: num
   :no-target:

******************
Performance Tuning
******************

This chapter discusses the various factors that impact GYRE's
execution time (i.e., how long it takes to complete a given
calculation), and presents strategies for minimizing this time.

Calculation Resolution
======================

GYRE's execution time depends on the resolution of the :ref:`spatial
<spatial-grids>` and :ref:`frequency <freq-grids>` grids it
uses. Specifically, for the :program:`gyre` frontend, the time to
process a single :nml:group:`mode` namelist group can be approximated
by

.. math::

   \tau \approx C_{\rm b} \, N \, \Nfreq + C_{\rm s} \, N \, \Nmode,

where :math:`N` is the number of spatial grid points, :math:`\Nfreq`
is the number of frequency grid points, :math:`\Nmode` is the number
of modes found, and :math:`C_{\rm b}` and :math:`C_{\rm s}` are
constants. The first (:math:`C_{\rm b}`) term represents the time
taken to bracket roots of the discriminant function, and the second
(:math:`C_{\rm s}`) the time taken to solve for these roots (see the
:ref:`numerical` chapter for details).

For the :program:`gyre_tides` frontend, the time to process a single
:nml:group:`orbit` namelist group can likewise be approximated by

.. math::

   \tau \approx C_{\rm t} \, N,

where :math:`C_{\rm t}` is a constant.

These expressions reveal that execution times can be shortened by
minimizing :math:`N` and :math:`\Nfreq` --- but this optimization
process must account for particular requirements for calculation
accuracy and completeness. :ads_citet:`townsend:2025` suggest an
strategy for adjusting :math:`N` to achieve a target accuracy, and the
:ref:`troubleshoot-miss` section explains how completeness can be
impacted by option choices in the :nml:group:`scan` namelist group.

.. nml:group:: scan
   :no-target:

A subtlety here involves the fact that :math:`N` is typically not
specified directly. When :nml:option:`scaffold_src
<grid.scaffold_src>` = :nml:value:`MODEL` (the default), the
baseline :math:`N` is established by the resolution of the stellar
model. However, the :ref:`iterative refinement algorithm
<spatial-grids-iter>` can increase :math:`N` by inserting additional
spatial grid points. This process implicitly depends on the
:nml:option:`freq_min` and :nml:option:`freq_max` options that
establish the range of the frequency scan; if this range includes
parts of the star's oscillation spectrum containing modes with very
large radial orders (whether p modes or g modes), then the refinement
algorithm will insert many points in order to resolve the modes'
wavefunctions. This can ultimately lead to huge :math:`N` and very
long execution times.

Calculation Numerics
====================

Another approach to reducing execution times is to change the
numerical algorithms employed by GYRE; in effect, this alters the
constantss :math:`C_{\rm b}`, :math:`C_{\rm s}` and :math:`C_{\rm t}`
appearing in the above expressions for :math:`\tau`.

Discretization Scheme
---------------------

.. nml:group:: num
   :no-target:

GYRE offers a number of schemes for :ref:`discretizing
<numerical-gyre-disc>` the oscillation equations. The choice of scheme
(selected using the :nml:option:`diff_scheme` option) represents a
trade-off between three factors: speed, accuracy and stability. In
general, higher-order discretization schemes are slower and
potentially less stable, but more accurate.

The performance of a given discretization scheme depends on many
factors, including the properties (e.g., CPU architecture, cache size,
memory bandwidth) of the hardware platform that GYRE runs
on. Therefore, empirical measurement remains the most practical
approach to characterizing behavior. :numref:`diff_scheme-table`
summarizes execution times for adiabatic test calculations\
[#test-diff_scheme]_ using the various choices for the
:nml:option:`diff_scheme` option. These data are calculated by
averaging over 10 separate runs on a single core of an `Intel Xeon
E5-2690 v4 CPU
<https://www.intel.com/content/www/us/en/products/sku/91770/intel-xeon-processor-e52690-v4-35m-cache-2-60-ghz/specifications.html>`__. Different
results will be obtained on other hardware, but the *relative*
performance of each scheme will likely remain similar.

.. _diff_scheme-table:

.. csv-table:: Execution time vs. diff_scheme
   :file: performance-tuning/stress.diff_scheme.csv
   :widths: auto
   :align: center

As discussed by :ads_citet:`townsend:2025`, the accuracy of
calculations involving :ref:`evolutionary <evol-models>` and
:ref:`polytropic <poly-models>` stellar models is typically limited by
the truncation error from interpolation, rather the truncation error
due to discretization. In such cases, there is little to be gained
from higher-order discretization schemes, and the
:nml:value:`COLLOC_GL2` scheme will remain the best choice for
performance-critical applications.

However, one notable exception to this recommendation concerns
non-adiabatic calculations. These are numerically challenging because
the oscillation equations deep in the stellar interior can be
extremely :wiki:`stiff <Stiff_equation>`. In such cases, the
:nml:value:`MAGNUS_GL2` scheme, although computationaly expensive,
appears to give superior results.

Matrix Solver
-------------

To solve the :ref:`linear system <numerical-gyre-linear>` formed by
the discretized oscillation equations, GYRE uses :wiki:`Gaussian
elimination` strategies that exploit the sparse matrix structure\
[#sparse-scaling]_. These solvers (selected using the
:nml:option:`ad_matrix_solver` and :nml:option:`nad_matrix_solver`
options) are numerically equivalent, but lead to different execution
times due to their differing approaches to accessing and manipulating
matrix elements.

:numref:`ad_matrix_solver-table` summarizes execution times for
adiabatic calculations\ [#test-matrix_solver]_ using the various
choices for the :nml:option:`ad_matrix_solver` option, on the same
Xeon CPU as before.

.. _ad_matrix_solver-table:

.. csv-table:: Execution times vs. ad_matrix_solver choice
   :file: performance-tuning/stress.ad_matrix_solver.csv
   :widths: auto
   :align: center

The winner is the :nml:value:`ROWPP` solver, recommending it as the
choice for all single-core calculations.

.. _performance-parallel:

Parallel Execution
==================

Parallel execution in GYRE is achieved using :wiki:`OpenMP`
multi-threading. For adiabatic calculations, GYRE parallelizes each
stage of the :ref:`mode search <numerical-gyre-search>`, yielding
speed-ups that scale with the thread count :math:`\Nthr`. This
is illustrated in :numref:`ad_matrix_solver-multi-table`, which
reprises :numref:`ad_matrix_solver-table` for different choices of
:math:`\Nthr`\ [#nthread]_. Note that the speed-up doesn't scale
exactly linearly with :math:`\Nthr`, due to :wiki:`dynamic frequency
scaling` of the CPU.

.. _ad_matrix_solver-multi-table:

.. csv-table:: Execution times vs. ad_matrix_solver choice for differing thread counts
   :file: performance-tuning/stress.ad_matrix_solver.multi.csv
   :widths: auto
   :align: center

For non-adiabatic calculations, the root-finding stage cannot be
parallelized when deflation is enabled (:nml:option:`deflate_roots` =
:nml:value:`.TRUE.`, the default), . The resulting performance loss
can partly be mitigated by using the :nml:value:`CYCLIC` matrix
solver, which is parallelized to yield a speed-up of order
:math:`\Nthr \log \Nthr`. This is illustrated in :numref:`table below
<matrix_solver-table>`, which reprises
:numref:`ad_matrix_solver-table` for non-adiabatic
calculations. Clearly, as :math:`\Nthr` increases, the
:nml:value:`CYCLIC` matrix solver overtakes :nml:value:`ROWPP` as the
fastest.

.. _nad_matrix_solver-multi-table:

.. csv-table:: Execution times vs. nad_matrix_solver choice for differing thread counts
   :file: performance-tuning/stress.nad_matrix_solver.multi.csv
   :widths: auto
   :align: center

To summarize this analysis:

* For single-threaded calculations, choose
  :nml:option:`ad_matrix_solver` = :nml:value:`ROWPP` and
  :nml:option:`nad_matrix_solver` = :nml:value:`ROWPP`.
* For multi-threaded calculations, choose
  :nml:option:`ad_matrix_solver` = :nml:value:`ROWPP` and
  :nml:option:`nad_matrix_solver` = :nml:value:`CYCLIC`.

Edge-cases can certainly arise where these heuristics don't apply;
then, as here, empirical measurement can guide the way.

.. rubric:: Footnotes


.. [#test-diff_scheme] See the
                       :file:`$GYRE_DIR/test/stress/stress.diff_scheme/run_stress.sh` script
                       for details of how the test calculations are run.

.. [#sparse-scaling] This is why the expressions for :math:`\tau`
                     scale with :math:`N`, rather than the usual
                     :math:`N^{3}` scaling that applies to Gaussian
                     elimination.

.. [#test-matrix_solver] See the
                        :file:`$GYRE_DIR/test/stress/stress.matrix_solver/run_stress.sh` script
                        for details of how the test calculations are run.

.. [#nthread] As set by the :envvar:`OMP_NUM_THREADS` environment variable.
