.. _comp-ptrope-solution:

Solution Method
===============

Specification
-------------

The structure of a composite polytrope is specified completely by

* a set of :math:`\nreg` polytropic indices :math:`n_{i}`
* a set of :math:`\nreg-1` boundary coordinates :math:`z_{i-1/2}`
* a set of :math:`\nreg` density jumps :math:`\Delta_{i-1/2} \equiv \ln [\rho_{i}(z_{i-1/2})/\rho_{i-1}(z_{i-1/2}]`

Although the normalizing densities :math:`\rho_{i,0}` have so far
been left unspecified, it's convenient to choose them as the density
at the beginning of their respective regions.

Solution
--------

The :ref:`structure equations <comp-ptrope-eqns>` may be solved as
an initial value problem. In the first region (:math:`i=1`) this IVP
involves integrating the Lane-Emden equation :eq:`lane-emden` from the
center :math:`z=0` to the first boundary :math:`z=z_{3/2}`, with the
initial conditions

.. math::

   \left.
   \begin{gathered}
   \theta_{i} = 1, \\
   \theta'_{i} = 0, \\
   B_{1} = 1, \\
   t_{1} = 1
   \end{gathered}
   \right\} \quad \text{at}\ z=0

(here, :math:`t_{i} \equiv \rho_{i,0}/\rho_{1,0}`).

The IVP in the intermediate regions (:math:`i = 2,\ldots,\nreg-1`)
involves integrating from :math:`z=z_{i-1/2}` to :math:`z=z_{i+1/2}`,
with initial conditions established from the preceding region via

.. math::

   \left.
   \begin{gathered}
   \theta_{i} = 1, \\
   \theta'_{i} = \frac{n_{i-1} + 1}{n_{i} + 1} \frac{\theta_{i-1}^{n_{i-1}}}{\theta_{i}^{n_{i}}} \frac{t_{i}}{t_{i-1}} \, \theta'_{i-1}, \\
   B_{i} = \frac{n_{i-1} + 1}{n_{i} + 1} \frac{\theta_{i}^{n_{i}+1}}{\theta_{i-1}^{n_{i-1}+1}} \frac{t_{i}^{2}}{t_{i-1}^{2}} \, B_{i-1}, \\
   \ln t_{i} = \ln t_{i-1} + n_{i-1} \ln \theta_{i-1} - n_{i} \ln \theta_{i} + \Delta_{i-1/2}.
   \end{gathered}
   \right\} \quad \text{at}\ z=z_{i-1/2}

The IVP in the final region (:math:`i=\nreg`) involves integrating
from :math:`z_{\nreg-1/2}` until :math:`\theta_{\nreg} = 0`. This
point defines the stellar surface, :math:`z=z_{\rm s}`. For some
choices of :math:`n_{i}`, :math:`z_{i-1/2}` and/or
:math:`\Delta_{i-1/2}`, the point :math:`\theta=0` can arise in an
earlier region :math:`i = \nreg_{\rm t} < \nreg`; in such cases, the
model specification must be truncated to :math:`\nreg_{\rm t}`
regions.
