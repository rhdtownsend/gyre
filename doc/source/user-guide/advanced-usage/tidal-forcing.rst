.. _tidal:

Tidal Forcing
=============

This section discusses how to evaluate the stellar response (fluid
displacements and perturbations) to tidal forcing, using the
:program:`gyre_tides` frontend. The response data can be used to
calculate the secular rates-of-change of orbital elements, or to
synthesize a light curve for a tidally distorted star.

Overview
--------

As discussed in the :ref:`osc-tidal` section, the tidal gravitational
potential (:eq:`e:tidal-pot`) of an orbiting companion can be
expressed as a superposition of partial potentials
:math:`\PhiTlmk`. For a given :nml_g:`tide` namelist group appearing
in the namelist input file, :program:`gyre_tides` solves for the
response of the star to each term in the superposition.

Truncating the Sums
-------------------

Although the sums appearing in eqn. (:eq:`e:tidal-pot`) are formally
infinite, the terms with large harmonic degree :math:`\ell` and/or
orbital harmonic :math:`k` typically produce a negligible
response. :program:`gyre_tides` offers a couple of approaches for
truncating the sums by dropping these terms. The simplest is to set
limits on the maximum values of the indices, through the
:nml_n:`l_max`, :nml_n:`k_min` and :nml_n:`k_max` parameters of the
:nml_g:`tide` namelist group (if desired, minimum values can also be
set using the corresponding :nml_n:`l_min` and :nml_n:`k_min`
parameters).

A slightly more sophisticated approach is to set these parameters to
large-ish values (say, a few hundred), and then also set one or both
of the :nml_n:`y_T_thresh_abs` and :nml_n:`y_T_thresh_rel`
parameters. These establish a threshold on the magnitude of

.. math::

   \yT \equiv \frac{\tPhiTlmk}{GM/R}

for a given tidal partial potential :math:`\tPhiTlmk` (see
eqn. :eq:`e:tidal-part-pot`) to be included in calculations; if
:math:`|\yT|` does not meet this threshold, it is ignored.

Optimizing Grids
----------------

During the :ref:`iterative refinement <spatial-grids-iter>` process
used in setting up spatial grids, the refinement criteria are
evaluated for every partial tide under consideration. If the
co-rotating forcing frequency

.. math::

   \sigmac \equiv k \Oorb - m \Orot

associated with a specific partial tidal potential is small compared
to the dynamical frequency of the star, many levels of refinement will
occur. While this is exactly what one wants in oscillation
calculations (because low-frequency modes have short spatial
wavelengths), it often isn't necessary in tidal calculations because
the response of a star to low-frequency forcing is the long-wavelength
equilibrium tide.

One way of preventing over-refinement due to low-frequency partial
potentials is to set the :nml_n:`omega_c_thresh` parameter in the
:nml_g:`tide` namelist group. This establishes a threshold on the
dimensionless frequency :math:`\omegac \equiv \sigmac \,
\sqrt{R^{3}/GM}`; partial potentials with :math:`|\omegac|` below this
threshold are treated as static tides (:math:`\omegac=0`), and are not
considered during the iterative refinement process.

An alternative approach is to avoid iterative refinement altogether,
instead obtaining the spatial grid from an external file (see the
:nml_v:`'FILE'` option of the :nml_n:`scaffold_src` parameter, in the
:ref:`grid-params` section). This is the most flexible approach, but
creating a grid that will adequately resolve the response to each
partial potential requires some care.

Output Files
------------

:program:`gyre_tides` writes response data to summary and detail
files. One detail file is created for each partial potential
evaluated, and the summary file collects together global data for all
partial potentials across all :nml_g:`tide` namelist groups. The
:nml_v:`id` output item can be used to determine which group a given
response belongs to.

The following Python code demonstrates how the summary data might be
used to evaluate the secular rates-of-change of orbital semi-major
axis, eccentricity, and argument of periastron, and the stellar
torque. The expression for :code:`e_dot` mirrors eqn. (23) of
:ads_citet:`sun:2023`, and for :code:`J_dot` eqn. (25) `ibid.`

.. literalinclude:: secular-rates.py
		    

