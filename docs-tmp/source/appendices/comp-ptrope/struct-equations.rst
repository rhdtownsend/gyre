.. _comp-ptrope-eqns:

Structure Equations
===================

Lane-Emden Equation
-------------------

In the :math:`i`'th region, a composite polytrope satisfies the
equation of hydrostatic equilibrium

.. math::

   -\frac{1}{\rho_{i}} \deriv{P_{i}}{r} = \deriv{\Phi_{i}}{r}

Substituting in the polytropic equation-of-state :eq:`poly-eos` yields

.. math::

   \frac{(n_{i}+1) P_{i,0}}{\rho_{i,0}^{1+1/n_{i}}} \deriv{}{r} \left( \rho_{i}^{1/n_{i}} \right) = - \deriv{\Phi_{i}}{r},

which can then be integrated with respect to :math:`r` to give

.. math::

   \frac{(n_{i}+1)P_{i,0}}{\Phi_{i,0} \, \rho_{i,0}} \left( \frac{\rho_{i}^{1/n_{i}}}{\rho_{i,0}^{1/n_{i}}} - 1 \right) = \left( 1 - \frac{\Phi_{i}}{\Phi_{i,0}} \right).

Here, the constants of integration have been chosen so that
:math:`\Phi_{i} = \Phi_{i,0}` when :math:`\rho_{i} =
\rho_{i,0}`. Rearranging, the density follows as

.. math::

   \rho_{i} = \rho_{i,0} \, \theta_{i}^{n_{i}},

where the polytropic dependent variable is introduced as

.. math::

   \theta_{i} = \left[ \frac{\Phi_{i,0} \, \rho_{i,0}}{(n_{i} + 1) P_{i,0}} \left( 1 - \frac{\Phi_{i}}{\Phi_{i,0}} \right)  + 1 \right].

With these expressions, Poisson's equation

.. math::

   \frac{1}{r^{2}} \deriv{}{r} \left( r^{2} \deriv{P_{i}}{r} \right) = 4 \pi G \rho_{i}

is recast as

.. math::

   \frac{1}{r^{2}} \deriv{}{r} \left( r^{2} \deriv{\theta_{i}}{r} \right) = - \frac{1}{A_{i}} \theta_{i}^{n_{i}},

where

.. math::

   A_{i} \equiv \frac{(n_{i} + 1) P_{i,0}}{4 \pi G \rho_{i,0}^{2}}.

A change of variables to the polytropic independent variable :math:`z
\equiv A_{1}^{-1/2} r` results in the dimensionless form

.. math::
   :label: lane-emden

   \frac{1}{z^{2}} \deriv{}{z} \left( z^{2} \deriv{\theta_{i}}{z} \right) = - B_{i} \theta_{i}^{n_{i}},

where :math:`B_{i} \equiv A_{1}/A_{i}`. This can be regarded as a
generalization of the usual :wiki:`Lane-Emden equation
<Lane%E2%80%93Emden_equation>` to composite polytropes.

Continuity Relations
--------------------

At the boundary between adjacent regions, the pressure and interior
mass must be continuous. If :math:`z_{i-1/2}` denotes the coordinate
of the boundary between the :math:`i-1` and :math:`i` regions, then
these continuity relations are expressed as

.. math::

   \left.
   \begin{gathered}
   B_{i} = \frac{n_{i-1} + 1}{n_{i} + 1} \frac{\theta_{i}^{n_{i}+1}}{\theta_{i-1}^{n_{i-1}+1}} \frac{\rho_{i,0}^{2}}{\rho_{i-1,0}^{2}} \, B_{i-1}, \\
   \theta'_{i} = \frac{n_{i-1} + 1}{n_{i} + 1} \frac{\theta_{i-1}^{n_{i-1}+1}}{\theta_{i}^{n_{i}+1}} \frac{\rho_{i,0}}{\rho_{i-1,0}} \, \theta'_{i-1},
   \end{gathered}
   \right\} \quad \text{at} \ z = z_{i-1/2}

respectively.
