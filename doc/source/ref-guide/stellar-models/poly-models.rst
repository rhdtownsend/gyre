.. _poly-models:

Polytropic Models
=================

Polytropic models represent stars in hydrostatic equilibrium that
follow the polytropic equation of state,

.. math::

   P = K \rho^{(n+1)/n}

for *piecewise-constant* :math:`K` and polytropic index :math:`n`. The
formalism behind these models is described in the :ref:`comp-ptrope`
section. The :program:`build_poly` executable can be used to create a
polytropic model and write it to a file in :ref:`POLY format
<poly-file-format>`.

To use a polytropic model with GYRE, set the :nml_n:`model_type`
parameter in the :nml_g:`model` namelist group to :nml_v:`'POLY'`,
and the :nml_n:`file` parameter to the name of the POLY file.

