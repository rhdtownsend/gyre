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

   z = \left[ \frac{K_{1} \rho_{\rm c}^{(1-n_{1})/n_{1}} (n_{1} + 1)}{4 \pi G} \right]^{-1/2} \, r,

where :math:`\rho_{\rm c}` is the central density of the star. The
dependent variable :math:`\theta_{i}` is related to the density, as
discussed below.

The structure of a composite polytrope is found by integrating the
composite Lane-Emden equation

.. math::

   \frac{1}{z^{2}} \deriv{}{z} \left( z^{2} \deriv{\theta_{i}}{z} \right) = - B_{i} \theta_{i}^{n_{i}}

in each region. Here, the independent variable is

.. math::

   z = \left[ \frac{K_{1} \rho_{\rm c}^{(1-n_{1})/n_{1}} (n_{1} + 1)}{4 \pi G} \right]^{-1/2} \, r,

where :math:`\rho_{\rm c}` is the central density of the star. The
dependent variable :math:`\theta_{i}` is related to the density, as
discussed in the :ref:`comp-ptrope-phys` section.

Initial & Boundary Conditions
-----------------------------

In the first region, the initial conditions are

.. math::

   \theta_{1}(0) = 1, \qquad \theta_{1}'(0) = 0

and :math:`B_{1} = 1`. In the subsequent regions, the initial conditions are

.. math::

   \theta_{i}(z_{i-1/2}) = 1, \qquad
   \theta'_{i}(z_{i-1/2}) = \frac{n_{i-1}+1}{n_{i}}
   \frac{\left[ \theta_{i}(z_{k-1/2}) \right]^{n_{i}+1}}{\left[ \theta_{i-1}(z_{i-1/2}) \right]^{n_{i-1}+1}}
   \frac{t_{i}}{t_{i-1}} \, \theta'_{i-1}(z_{i-1/2})

where :math:`z_{i-1/2}` is the coordinate of the boundary between the
:math:`i-1` and :math:`i` regions; moreover,

.. math::

   B_{i} = B_{i-1} \frac{n_{i-1}+1}{n_{i}+1}
   \frac{\left[ \theta_{i}(z_{k-1/2}) \right]^{n_{i}+1}}{\left[ \theta_{i-1}(z_{i-1/2}) \right]^{n_{i-1}+1}}
   \frac{t_{i}^{2}}{t_{i-1}^{2}}.

The recurrence for :math:`t_{i}` is

.. math::

   \ln t_{k} = \ln t_{i-1} + n_{i-1} \ln \theta_{i-1}(z_{i-1/2}) - n_{i} \ln \theta_{i}(z_{i-1/2}) + \Delta_{i-1/2},

with :math:`t_{1} = 1` and

.. math::

   \Delta_{i-1/2} = \ln \left[ \frac{\rho_{i}(z_{i-1/2})}{\rho_{i-1}(z_{i-1/2})} \right]

quantifying the density jump at :math:`z_{i-1/2}`.

The surface of the composite polytropic model, :math:`z=z_{\rm
s}`, is defined implicitly by the boundary condition

.. math::

   \theta_{N}(z_{\rm s}) = 0.
