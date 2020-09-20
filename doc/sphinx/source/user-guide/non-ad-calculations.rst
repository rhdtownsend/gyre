.. _non-ad-calcs:

**************************
Non-Adiabatic Calculations
**************************

Overview
========

In addition to calculating the frequencies/periods of oscillation
modes, GYRE can obtain the corresponding linear growth or damping
rates. To do this, it includes *non-adiabatic* terms in the
:ref:`oscillation equations <osc-equations>` describing the transfer of heat
between neighboring fluid elements. With these terms, the equations
and their solutions become complex quantities. The assumed time
dependence for perturbations is :math:`\propto \exp (-\ii \sigma t)`
(see the :ref:`linear-equations` section); therefore, the real part
:math:`\sigmar` and imaginary part :math:`\sigmai` of the
eigenfrequency are related to the mode period :math:`\Pi` and growth
e-folding time :math:`\tau`, respectively, via

.. math::

   \Pi = \frac{2\pi}{\sigmar}, \qquad
   \tau = \frac{1}{\sigmai}.

Non-adiabatic calculations proceed using the same general approach as
in the adiabatic case, by searching for the roots of a discriminant
function :math:`\Dfunc(\omega)` (see the :ref:`gyre-fundamentals`
chapter for more details). However, a challenge is that there is no
simple way to bracket roots in the complex plane. Instead, GYRE must
generate initial trial roots that are close to the true roots, and
then refine them iteratively. Currently, GYRE offers three methods for
establishing the trial roots.

Adiabatic Method
================

The adiabatic method involves adopting the (real) roots found from
adiabatic calculations as the initial trial roots for the
non-adiabatic problem. This works well as long as the two sets of
roots lie close together in the complex plane --- typically, when the
oscillation modes are only weakly non-adiabatic, with
:math:`|\sigmai/\sigmar| \ll 1`.

To perform non-adiabatic calculations with the adiabatic method, set
the following parameters in the :nml_g:`osc` namelist group:

* :nml_n:`adiabatic`\ =\ :nml_v:`.TRUE.`\ [#default]_
* :nml_n:`nonadiabatic`\ =\ :nml_v:`.TRUE.`

and the following parameters in the :nml_g:`num` namelist group:

* :nml_n:`ad_search`\ =\ :nml_v:`'BRACKET'`\ [#default]_
* :nml_n:`nad_search`\ =\ :nml_v:`'AD'`

You may also wish to use the following setting in the :nml_g:`num`
namelist group:

* :nml_n:`diff_scheme`\ =\ :nml_v:`'MAGNUS_GL2'`

This tells GYRE to evaluate the finite-difference equations using the
2nd order Magnus scheme; experience suggests that this gives the most
reliable convergence for the root refinement.

Minmod Method
=============

The minmod method involves evaluating the discriminant function along
the real-:math:`\omega` axis, and then adopting local minima in its
modulus :math:`|\Dfunc|` as the initial trial roots for the
non-adiabatic problem. The method is described in full in
:ads_citet:`goldstein:2020`; as shown there, it does not perform
significantly better than the adiabatic method, and is included in
GYRE for the sake of completeness.

To perform non-adiabatic calculations with the adiabatic method, set
the following parameters in the :nml_g:`osc` namelist group:

* :nml_n:`adiabatic`\ =\ :nml_v:`.FALSE.`\ [#optional]_
* :nml_n:`nonadiabatic`\ =\ :nml_v:`.TRUE.`

and the following parameters in the :nml_g:`num` namelist group:

* :nml_n:`nad_search`\ =\ :nml_v:`'MINMOD'`

As with the adiabatic method, you may also wish to use the following
setting in the :nml_g:`num` namelist group:

* :nml_n:`diff_scheme`\ =\ :nml_v:`'MAGNUS_GL2'`

Contour Method
==============

The contour method involves evaluating the discriminant function on a
grid in the complex-:math:`\omega` plane, and then adopting
intersections between the real zero-contours :math:`\Dfuncr=0`, and
the corresponding imaginary ones :math:`\Dfunci=0`, as the initial
trial roots for the non-adiabatic problem. The method is described in
full in :ads_citet:`goldstein:2020`; it is very effective even for
strongly non-adiabatic modes with :math:`|\sigmai/\sigmar| \sim 1`,
although there is an increased computational cost (see :ref:`here <faq-cluster>`
for one strategy for mitigating this cost).

To perform non-adiabatic calculations with the contour method, set
the following parameters in the :nml_g:`osc` namelist group:

* :nml_n:`adiabatic`\ =\ :nml_v:`.FALSE.`\ [#optional]_
* :nml_n:`nonadiabatic`\ =\ :nml_v:`.TRUE.`

and the following parameters in the :nml_g:`num` namelist group:

* :nml_n:`nad_search`\ =\ :nml_v:`'CONTOUR'`

Finally, you must also ensure that at least one :nml_g:`scan` namelist
group with :nml_n:`axis`\ =\ :nml_v:`'REAL'` is present, and likewise
at least one with :nml_n:`axis`\ =\ :nml_v:`'IMAG'`. Together, these
groups define the real and imaginary axes of the discriminant grid in
the complex-:math:`\omega` plane. As a rule of thumb, the resolution
along the imaginary axis should be comparable to that along the real
axis; this ensures that the contour-tracing algorithm behaves well.

Finally, as with the adiabbatic method, you may also wish to use the
following setting in the :nml_g:`num` namelist group:

* :nml_n:`diff_scheme`\ =\ :nml_v:`'MAGNUS_GL2'`

.. rubric:: Footnotes

.. [#default] This is the default setting; you don't need to include it explicitly

.. [#optional] This is optional; leave it out if you want GYRE to perform adiabatic calculations as well
