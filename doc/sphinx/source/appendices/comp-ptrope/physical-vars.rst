.. _comp-ptrope-phys:

Physical Variables
==================

Once the :ref:`Lane-Emden equation <comp-ptrope-lane-emden>` has been solved, the density in each
region can be evaluated by

.. math::

   \rho_{i} = \rho_{\rm c} t_{i} \theta_{i}^{n_{i}}.

The pressure then follows from the :ref:`polytropic equation of state <comp-ptrope-eos>`.

The interior mass :math:`m` is evaluated by introducing the auxillary
quantity :math:`\mu`, which is defined in the first region by

.. math::

   \mu_{1}(z) = - z^{2} \theta'_{1} (z),

and in subsequent regions by

.. math::

   \mu_{i}(z) = \mu_{i-1}(z_{i-1/2}) - \frac{t_{i}}{B_{i}} \left[ z^{2} \theta'_{i} (z) - z_{i-1/2}^{2} \theta'_{i} (z_{i-1/2}) \right].

The interior mass then follows as

.. math::

   m = M \frac{\mu_{i}(z)}{\mu_{\rm s}},

where :math:`\mu_{\rm s} \equiv \mu_{N}(z_{\rm s})`.

