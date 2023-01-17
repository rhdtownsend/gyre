.. _tidal-sep-eqns:

Separated Equations
===================

Because the tidal potential (:eq:`e:tidal-pot`) superposes many
different spherical harmonics, the separation of variables
(:eq:`e:osc-sol-forms`) applied to the oscillation equations must
be replaced by the more-general expressions

.. math::
   :label: e:tidal-sol-forms

   \begin{aligned}
   \xir(r,\theta,\phi;t) &= \sum_{\ell,m,k} \txirlmk(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii k \Oorb t), \\
   \xit(r,\theta,\phi;t) &= \sum_{\ell,m,k} \txihlmk(r) \, \pderiv{}{\theta} Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \Oorb t), \\
   \xip(r,\theta,\phi;t) &= \sum_{\ell,m,k} \txihlmk(r) \, \frac{\ii m}{\sin\theta} Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii k \Oorb t), \\
   f'(r,\theta,\phi;t) &= \sum_{\ell,m,k} \tflmk'(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii k \Oorb t)
   \end{aligned}

(the notation for the sums has been abbreviated). Substituting these
solution forms into the :ref:`linearized equations
<tidal-linear-eqns>`, and taking advantage of the orthonormality of
the spherical harmonics, leads to a set of differential equations for
each combination of :math:`\ell`, :math:`m` and :math:`k`. A given set
resembles the corresponding :ref:`oscillation equations
<osc-sep-eqns>`, with just a couple changes:

- Rather than being an eigenvalue parameter, the oscillation frequency
  is set by :math:`\sigma = k \Oorb`, representing the forcing
  frequency of the partial tidal potential in an inertial frame.
- The perturbation :math:`\tPhi'` is replaced by :math:`\tPsi' \equiv
  \tPhi' + \tPhiT`, representing the radial part of the total (self +
  tidal) gravitational potential perturbation.
