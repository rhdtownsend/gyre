.. _comp-ptrope-vars:

Physical Variables
==================

Once the Lane-Emden equation :eq:`lane-emden` has been solved, the density in each
region can be evaluated by

.. math::

   \rho_{i} = \rho_{1,0} \, t_{i} \, \theta_{i}^{n_{i}}.

The pressure then follows from the equation-of-state
:eq:`poly-eos` as

.. math::

   P_{i} = P_{1,0} \, \frac{n_{1}+1}{n_{i}+1} \, \frac{t_{i}^{2}}{B_{i}} \, \theta_{i}^{n_{i}+1}.

The interior mass :math:`m` is evaluated by introducing the auxillary
quantity :math:`\mu`, which is defined in the first region by

.. math::

   \mu_{1}(z) = - z^{2} \theta'_{1} (z),

and in subsequent regions by

.. math::

   \mu_{i}(z) = \mu_{i-1}(z_{i-1/2}) - \frac{t_{i}}{B_{i}} \left[ z^{2} \theta'_{i} (z) - z_{i-1/2}^{2} \theta'_{i} (z_{i-1/2}) \right].

The interior mass then follows as

.. math::

   M_{r} = M \frac{\mu_{i}}{\mu_{\rm s}},

where :math:`\mu_{\rm s} \equiv \mu_{\nreg}(z_{\rm s})`.

