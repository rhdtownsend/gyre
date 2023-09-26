.. _osc-sep-eqns:

Separated Equations
===================

With a separation of variables in spherical-polar coordinates
:math:`(r,\theta,\phi)`, and assuming an oscillatory time (:math:`t`)
dependence with angular frequency :math:`\sigma`, solutions to the
:ref:`linearized equations <osc-linear-eqns>` can be expressed as

.. math::
   :label: e:osc-sol-forms

   \begin{aligned}
   \xir(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \txir(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t) \right], \\
   \xit(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \txih(r) \, \pderiv{}{\theta} Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t) \right], \\
   \xip(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \txih(r) \, \frac{\ii m}{\sin\theta} Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t) \right], \\
   f'(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \tf'(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t) \right].
   \end{aligned}

Here, :math:`\xir`, :math:`\xit` and :math:`\xip` are the radial,
polar and azimuthal components of the displacement perturbation vector
:math:`\vxi`; :math:`Y^{m}_{\ell}` is the spherical harmonic with
harmonic degree :math:`\ell` and azimuthal order :math:`m`; and again
:math:`f` stands for any perturbable scalar. The displacement
perturbation vector is related to the velocity perturbation via

.. math::

   \vv' = \pderiv{\vxi}{t}.

Substituting the above solution forms into the linearized equations,
the mechanical (mass and momentum conservation) equations become

.. math::
   :label: e:osc-sep-cont

   \trho' + \frac{1}{r^{2}} \deriv{}{r} \left( \rho r^{2} \txir \right) - \frac{\ell(\ell+1)}{r} \rho \txih = 0,

.. math::

   -\sigma^{2} \rho \txir = - \deriv{\tP'}{r} - \trho' \deriv{\Phi}{r} - \rho \deriv{\tPhi'}{r},

.. math::

   -\sigma^{2} \rho r \txih = - \tP' - \rho \tPhi'.

Likewise, Poisson's equation becomes

.. math::
   :label: e:osc-sep-poisson

   \frac{1}{r^{2}} \deriv{}{r} \left( r^{2} \deriv{\tPhi'}{r} \right) - \frac{\ell(\ell+1)}{r^{2}} \txih = 4 \pi G \trho'

and the heat equation becomes

.. math::

   -\ii \sigma T \delta \tS = \delta \tepsnuc
   - \deriv{\delta \tLrad}{M_{r}} + \frac{\ell(\ell+1)}{\sderiv{\ln T}{r}} \frac{\Fradr}{\rho}  \frac{\tT'}{T} +
    \ell(\ell + 1) \frac{\txih}{r} \deriv{\Lrad}{M_{r}},

where

.. math::

   \delta \tLrad \equiv 4 \pi r^{2} \left( \delta \tFradr + 2 \frac{\txir}{r} \Fradr \right)

is the Lagrangian perturbation to the radiative luminosity. The radial part of the radiative diffusion equation becomes

.. math::

   \delta\tFradr = \left[
   4 \frac{\delta \tT}{T} - \frac{\delta\trho}{\rho} - \frac{\delta\tkappa}{\kappa} + 
   \frac{\sderiv{(\delta \tT/T)}{\ln r}}{\sderiv{\ln T}{\ln r}} \right] \Fradr.

Finally, the thermodynamic, nuclear and opacity relations become

.. math::

   \frac{\delta \trho}{\rho} = \frac{1}{\Gammi} \frac{\delta \tP}{P} - \upsT \frac{\delta \tS}{\cP},
   \qquad
   \frac{\delta \tT}{T} = \nabla_{\rm ad} \frac{\delta \tP}{P} + \frac{\delta \tS}{\cP},

.. math::

   \frac{\delta \tepsnuc}{\epsnuc} = \epsnucad \frac{\delta \tP}{P} + \epsnucS \frac{\delta \tS}{\cP},
   \qquad
   \frac{\delta \tkappa}{\kappa} = \kapad \frac{\delta \tP}{P} + \kapS \frac{\delta \tS}{\cP}.

