.. _comp-ptrope-le:

Lane-Emden Equation
===================

Differential Equation
---------------------

The structure of a composite polytrope is found by integrating the
composite Lane-Emden equation

.. math::

   \frac{1}{z^{2}} \deriv{}{z} \left( z^{2} \deriv{\theta_{i}}{z} \right) = - B_{i} \theta_{i}^{n_{i}}

in each region. Here, the independent variable is

.. math::

   z = \left[ \frac{K_{1} (\rhoc t_{1})^{(1-n_{1})/n_{1}} (n_{1} + 1)}{4 \pi G} \right]^{-1/2} \, r,

where :math:`\rhoc` is the central density of the star; likewise, the
dependent variable is

.. math::

   \theta_{i} = \left( \frac{1}{t_{i}} \frac{\rho_{i}}{\rhoc} \right)^{1/n_{i}},

where :math:`t_{i}` is a term introduced to allow for density
discotinuities at the boundaries between regions. The quantity
:math:`B_{i}` depends on :math:`K_{1}` and :math:`K_{i}` via

.. math::

   B_{i} = \frac{K_{1}}{K_{i}} \frac{n_{1}+1}{n_{i}+1} \frac{(\rhoc t_{1})^{(1-n_{1})/n_{1}}}{(\rhoc t_{i})^{(1-n_{i})/n_{i}}};

clearly, :math:`B_{1} = 1`.
 
Initial & Boundary Conditions
-----------------------------

In the first region, the initial conditions are

.. math::

   \theta_{1}(0) = 1, \qquad \theta_{1}'(0) = 0.

In the subsequent regions, the initial conditions are

.. math::

   \theta_{i}(z_{i-1/2}) = 1, \qquad
   \theta'_{i}(z_{i-1/2}) = \frac{n_{i-1}+1}{n_{i}}
   \frac{\left[ \theta_{i}(z_{k-1/2}) \right]^{n_{i}+1}}{\left[ \theta_{i-1}(z_{i-1/2}) \right]^{n_{i-1}+1}}
   \frac{t_{i}}{t_{i-1}} \, \theta'_{i-1}(z_{i-1/2})

where :math:`z_{i-1/2}` is the coordinate of the boundary between the
:math:`i-1` and :math:`i` regions. Pressure continuity at this
boundary requires that

.. math::

   B_{i} = B_{i-1} \frac{n_{i-1}+1}{n_{i}+1}
   \frac{\left[ \theta_{i}(z_{k-1/2}) \right]^{n_{i}+1}}{\left[ \theta_{i-1}(z_{i-1/2}) \right]^{n_{i-1}+1}}
   \frac{t_{i}^{2}}{t_{i-1}^{2}}.

The recurrence for :math:`t_{i}` is

.. math::

   \ln t_{k} = \ln t_{i-1} + n_{i-1} \ln \theta_{i-1}(z_{i-1/2}) - n_{i} \ln \theta_{i}(z_{i-1/2}) + \Delta_{i-1/2},

where :math:`t_{1} = 1`, and

.. math::

   \Delta_{i-1/2} = \ln \left[ \frac{\rho_{i}(z_{i-1/2})}{\rho_{i-1}(z_{i-1/2})} \right]

quantifies the density jump at :math:`z_{i-1/2}`.

The surface of the composite polytropic model, :math:`z=z_{\rm
s}`, is defined implicitly by the boundary condition

.. math::

   \theta_{\nreg}(z_{\rm s}) = 0.
