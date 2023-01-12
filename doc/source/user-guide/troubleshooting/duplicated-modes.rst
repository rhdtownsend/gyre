.. _troubleshoot-dupe:

Duplicated Modes
================

Sometimes two oscillation modes with the same :math:`\npg` are found
during a :program:`gyre` calculation. This violates the expectation
that :math:`\npg` be monotonic-increasing, and can happen for a few
reasons.

Bad Stellar Model
-----------------

If the input stellar doesn't conserve mass properly, then one or more
bogus (unphysical) modes can appear with the same :math:`\npg` as an
existing mode, but a significantly different frequency. Such modes can
also be spotted because their radial order is very different from the
adjacent-in-frequency modes.

Bogus modes arise because the input stellar model doesn't conserve
mass. GYRE` assumes that the density :math:`\rho` and interior mass
:math:`M_{r}` are related by equation :eq:`mass-eq`. Given that there
are many different ways to discretize this equation, there is a
certain amount of numerical 'slop' that arises when going from
discrete form assumed in the stellar evolution code that generated the
model (e.g., MESA) through to the discrete form implicitly assumed by
GYRE. This slop can sometimes produce bogus modes, and the fix is to
re-create the model with a higher spatial resolution.

Another possibility arises if the input stellar model contains density
discontinuities that aren't properly marked using :ref:`double points
<evol-models-double>`. Often, these discontinuities can be diagnosed
by plotting the :math:`I_0` or :math:`I_1` first integrals (for radial
and dipole modes, respectively) as a function of radius. These first
integrals should vanish everywhere (as shown by
:ads_citealt:`takata:2006a`), but will typically show abrupt jumps to
non-zero values at the location of unmarked discontinuities. The fix
is to re-create the model with double points inserted as necessary; in
the case of MESA models, this can be achieved using the
:nml_n:`add_double_points_to_pulse_data` parameter of the
:nml_g:`controls` namelist group.

Non-adiabatic Effects
---------------------

When non-adiabatic effects cause a mode to be mis-classified (as
discussed in the :ref:`troubleshoot-miss` section), often the incorrect
:math:`\npg` value duplicates that of another mode. As before, a
mis-classification must be fixed manually by determining which
adiabatic mode the problematic non-adiabatic mode corresponds to.
