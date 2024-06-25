.. _non-ad-osc:

Non-Adiabatic Oscillations
==========================

This section discusses how to undertake *non-adiabatic* oscillation
calculations using the :program:`gyre` frontend. Asteroseismic studies
typically rely on adiabatic calculations, because the frequencies of
oscillation modes are the primary focus. However, for heat-driven
modes the linear growth or damping rates can also be of interest ---
and evaluating these requires that non-adiabatic effects are included
in the oscillation equations.

.. note::
   Not all types of stellar mode include the necessary data
   (e.g., thermodynamic coefficients, opacity partial derivatives) to
   undertake non-adiabatic calculations. The :ref:`model-caps` section
   summarizes this information.

Overview
--------

To include non-adiabatic effects :program:`gyre` augments the
linearized mass, momentum and Poisson equations with the linearized
heat and radiative diffusion equations (see the :ref:`osc-linear-eqns`
section for full details). With these additions, the equations and
their solutions become complex quantities. The assumed time dependence
for perturbations is :math:`\propto \exp (-\ii \sigma t)`; therefore,
the real part :math:`\sigmar` and imaginary part :math:`\sigmai` of
the eigenfrequency are related to the mode period :math:`\Pi` and
growth e-folding time :math:`\tau`, respectively, via

.. math::

   \Pi = \frac{2\pi}{\sigmar}, \qquad
   \tau = \frac{1}{\sigmai}.

Solving the non-adiabatic equations proceeds using the same general
approach as in the adiabatic case, by searching for the roots of a
discriminant function :math:`\Dfunc(\omega)` (see the :ref:`numerical`
chapter for more details). However, a challenge is that there is no
simple way to bracket roots in the complex plane. Instead,
:program:`gyre` must generate initial trial roots that are close to
the true roots, and then refine them iteratively. Currently,
:program:`gyre` offers three methods for establishing the trial roots.

.. _non-ad-ad:

Adiabatic Method
----------------

The adiabatic method involves adopting the (real) roots found from
adiabatic calculations as the initial trial roots for the
non-adiabatic problem. This works well as long as the adiabatic and
non-adiabatic roots lie close together in the complex plane ---
typically, when the oscillation modes are only weakly non-adiabatic,
with :math:`|\sigmai/\sigmar| \ll 1`.

To perform non-adiabatic calculations with the adiabatic method, set
the following parameters in the :nml_g:`osc` namelist group:

* :nml_n:`nonadiabatic`\ =\ :nml_v:`.TRUE.`
* :nml_n:`adiabatic`\ =\ :nml_v:`.TRUE.`\ [#default]_

and the following parameters in the :nml_g:`num` namelist group:

* :nml_n:`ad_search`\ =\ :nml_v:`'BRACKET'`\ [#default]_
* :nml_n:`nad_search`\ =\ :nml_v:`'AD'`\ [#default]_

You may also wish to use the following setting in the :nml_g:`num`
namelist group:

* :nml_n:`diff_scheme`\ =\ :nml_v:`'MAGNUS_GL2'`

This tells :program:`gyre` to evaluate the finite-difference equations
using the 2nd order Magnus scheme; experience suggests that this gives
the most reliable convergence for the root refinement.

An example of the adiabatic method in action can be found in the
:file:`{$GYRE_DIR}/test/nad/mesa/bcep/gyre.in` namelist input file,
which is set up to find :math:`\ell=0,\ldots,3` modes of a
:math:`20\,\Msun` :math:`\beta` Cephei model using the adiabatic
method. The important parts are as follows:

.. literalinclude:: non-ad-oscillations/gyre.in
   :language: console
   :start-at: &osc
   :end-before: &grid

Note the :nml_n:`nonadiabatic` parameter in the :nml_g:`osc` namelist
group, and the :nml_n:`diff_scheme` parameter in the :nml_g:`num`
namelist group. The :nml_n:`restrict_roots`\ =\ :nml_v:`.FALSE.`
setting, also in the :nml_g:`num` namelist group, tells
:program:`gyre` not to reject any modes that have :math:`\sigmar`
outside the frequency range specified by the :nml_g:`scan` namelist
group; this ensures that modes whose non-adiabatic frequencies fall
just outside the frequency grid are still found.

.. _non-ad-minmod:

Minmod Method
-------------

The minmod method involves evaluating the discriminant function along
the real-:math:`\omega` axis, and then adopting local minima in its
modulus :math:`|\Dfunc|` as the initial trial roots for the
non-adiabatic problem. The method is described in full in
:ads_citet:`goldstein:2020`; as shown there, it does not perform
significantly better than the adiabatic method, and is included in
:program:`gyre` for the sake of completeness.

To perform non-adiabatic calculations with the minmod method, set
the following parameters in the :nml_g:`osc` namelist group:

* :nml_n:`nonadiabatic`\ =\ :nml_v:`.TRUE.`
* :nml_n:`adiabatic`\ =\ :nml_v:`.FALSE.`\ [#optional]_

and the following parameters in the :nml_g:`num` namelist group:

* :nml_n:`nad_search`\ =\ :nml_v:`'MINMOD'`

As with the adiabatic method, you may also wish to use the following
setting in the :nml_g:`num` namelist group:

* :nml_n:`diff_scheme`\ =\ :nml_v:`'MAGNUS_GL2'`

An example of the minmod method in action can be found in the
:file:`{$GYRE_DIR}/test/nad/mesa/bcep-minmod/gyre.in` namelist input
file, which is equivalent to
:file:`{$GYRE_DIR}/test/nad/mesa/bcep/gyre.in` but using the
minmod method. The important parts are as follows:

.. literalinclude:: non-ad-oscillations/gyre-minmod.in
   :language: console
   :start-at: &osc
   :end-before: &grid

Note the additional :nml_n:`nad_search`\ =\ :nml_v:`'MINMOD'` parameter
in the :nml_g:`num` namelist group, which stipulates that the minmod
method should be used.

.. _non-ad-contour:

Contour Method
--------------

The contour method involves evaluating the discriminant function on a
grid in the complex-:math:`\omega` plane, and then adopting
intersections between the real zero-contours :math:`\Dfuncr=0`, and
the corresponding imaginary ones :math:`\Dfunci=0`, as the initial
trial roots for the non-adiabatic problem. The method is described in
full in :ads_citet:`goldstein:2020`; it is very effective even for
strongly non-adiabatic modes with :math:`|\sigmai/\sigmar| \sim 1`,
although there is an increased computational cost.

To perform non-adiabatic calculations with the contour method, set
the following parameters in the :nml_g:`osc` namelist group:

* :nml_n:`nonadiabatic`\ =\ :nml_v:`.TRUE.`
* :nml_n:`adiabatic`\ =\ :nml_v:`.FALSE.`\ [#optional]_

and the following parameters in the :nml_g:`num` namelist group:

* :nml_n:`nad_search`\ =\ :nml_v:`'CONTOUR'`

You must also ensure that at least one :nml_g:`scan` namelist
group with :nml_n:`axis`\ =\ :nml_v:`'REAL'` is present, and likewise
at least one with :nml_n:`axis`\ =\ :nml_v:`'IMAG'`. Together, these
groups define the real and imaginary axes of the discriminant grid in
the complex-:math:`\omega` plane. As a rule of thumb, the resolution
along the imaginary axis should be comparable to that along the real
axis; this ensures that the contour-tracing algorithm behaves well.

Finally, as with the adiabatic method, you may also wish to use the
following setting in the :nml_g:`num` namelist group:

* :nml_n:`diff_scheme`\ =\ :nml_v:`'MAGNUS_GL2'`

.. note::

   Because g modes are spaced uniformly in period (in the asymptotic
   limit of large radial order), it would seem sensible to set
   :nml_n:`grid_type`\ =\ :nml_v:`'INVERSE'` in the :nml_g:`scan`
   namelist group(s) that correspond to the real axis (i.e.,
   :nml_n:`axis`\ =\ :nml_v:`'REAL'`). However, this typically results
   in a mismatch between the resolution of the real and imaginary
   axes, and the contour method doesn't perform well. A fix for this
   issue will be forthcoming in a future release of GYRE,
   but in the meantime it's probably best to avoid the contour method
   for g modes.

An example of the minmod method in action can be found in the
:file:`{$GYRE_DIR}/test/nad/mesa/bcep-contour/gyre.in` namelist input
file, which is equivalent to
:file:`{$GYRE_DIR}/test/nad/mesa/bcep/gyre.in` but using the
minmod method. The important parts are as follows:

.. literalinclude:: non-ad-oscillations/gyre-contour.in
   :language: console
   :start-at: &osc
   :end-before: &grid

Note the additional :nml_n:`nad_search`\ =\ :nml_v:`'CONTOUR'`
parameter in the :nml_g:`num` namelist group, which stipulates that
the contour method should be used; and, the fact that there are now
two :nml_g:`scan` namelist groups, one with :nml_n:`axis`\ =\
:nml_v:`'REAL'` and the other with :nml_n:`axis`\ =\ :nml_v:`'IMAG'`.

.. rubric:: Footnotes

.. [#default] This is the default setting; you don't need to include it explicitly.

.. [#optional] This is optional; leave it out if you want :program:`gyre` to perform adiabatic calculations as well.
