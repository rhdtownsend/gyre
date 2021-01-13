.. _comp-ptrope-eos:

Equation of State
=================

Consider a composite polytrope composed of :math:`\nreg` regions
extending from the origin out to the stellar surface. In the
:math:`i`'th region (:math:`1 \leq i \leq \nreg`), the pressure
:math:`P` and density :math:`\rho` are related by the polytropic
equation-of-state

.. math::
   :label: poly-eos

   \frac{P_{i}}{P_{i,0}} = \left( \frac{\rho_{i}}{\rho_{i,0}} \right)^{(n_{i} + 1)/n_{i}}

where the normalizing pressure :math:`P_{i,0}` and density
:math:`\rho_{i,0}`, together with the polytropic index :math:`n_{i}`,
are constant across the region but may change from one region to the
next. At the :math:`\nreg-1` boundaries between adjacent regions, the
pressure and interior mass :math:`M_{r}` are required to be
continuous, but the density may jump.
