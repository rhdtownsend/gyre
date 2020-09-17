.. _poly-models:

Polytropic Models
=================

Polytropic models represent stars in hydrostatic equilibrium that
follow the polytropic equation of state,

.. math::

   p = K \rho^{(\npoly+1)/\npoly},

for *piecewise-constant* :math:`K` and polytropic index :math:`\npoly`. The structure of
polytropic models is obtained by solving the Lane-Emden equation; GYRE
doesn't do this itself, but can read structure data from external
files created by the :program:`build_poly` executable (see the
:ref:`build-poly` appendix for more details).

To use a polytropic model with GYRE, set the :nml_n:`model_type`
parameter in the :nml_g:`model` namelist group to :nml_v:`'POLY'`,
and the :nml_n:`file` parameter to the name of the polytropic
structure file.

