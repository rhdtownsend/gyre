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
then refine them iteratively. Currently, GYRE offers two alternative
methods for establishing the trial roots.

Adiabatic Method
================

The adiabatic method involves adopting the (real) roots found from
adiabatic calculations as the initial trial roots for the
non-adiabatic problem. This works well as long as the two sets of
roots lie close together in the complex plane --- typically, when the
oscillation modes are only weakly non-adiabatic, with
:math:`|\sigmai/\sigmar| \ll 1`.

To perform non-adiabatic calculations with the adiabatic method, use
the standard :program:`gyre` executable with the following parameter
settings in the :nml_g:`osc` namelist group:

* :nml_n:`adiabatic`\ =\ :nml_v:`.TRUE.`
* :nml_n:`nonadiabatic`\ =\ :nml_v:`.TRUE.`

You may also wish to use the following setting in the :nml_g:`num`
namelist group:

* :nml_n:`diff_scheme`\ =\ :nml_v:`'MAGNUS_GL2'`

This tells GYRE to evaluate the finite-difference equations using the
2nd order Magnus scheme; experience suggests that this gives the most
reliable convergence for the root refinement.
  
Contour Method
==============

The contour method involves evaluating the discriminant function on a
grid in the complex-:math:`\omega` plane, and then adopting
intersections between the real zero-contours :math:`\Dfuncr=0`, and
the corresponding imaginary ones :math:`\Dfunci=0`, as the initial
trial roots for the non-adiabatic problem. The method is described in
full in :ads_citet:`goldstein:2020`; while it is computationally
expensive, it performs well even for strongly non-adiabatic modes with
:math:`|\sigmai/\sigmar| \sim 1`.

To perform non-adiabatic calculations with the contour method, use the
:program:`gyre_contour` executable. Broadly, this works in the same
way as the the :program:`gyre` executable: run parameters are read
from an input file defining a number of Fortran namelist
groups. However, an important difference is that at least one
:nml_g:`scan` namelist group with :nml_n:`axis`\ =\ :nml_v:`'REAL'`
must be present, and likewise at least one with :nml_n:`axis`\ =\
:nml_v:`'IMAG'`. Together, these groups define the real and imaginary
axes of the discriminant grid in the complex-:math:`\omega` plane.

As a rule of thumb, the resolution along the imaginary axis should be
comparable to that along the real axis; this ensures that the
contour-tracing algorithm behaves well. You may also wish to use the
following setting in the :nml_g:`num` namelist group:

* :nml_n:`diff_scheme`\ =\ :nml_v:`'MAGNUS_GL2'`

This tells GYRE to evaluate the finite-difference equations using the
2nd order Magnus scheme; experience suggests that this gives the most
reliable convergence for the root refinement.
