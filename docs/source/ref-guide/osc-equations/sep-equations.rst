.. _osc-sep-eqns:

Separated Equations
===================

The :ref:`linearized equations <osc-linear-eqns>` are fully separable
in spherical-polar coordinates :math:`(r,\theta,\phi)` and time
:math:`t`. Assuming that all perturbations have an oscillatory time
dependence :math:`\propto \exp(-\ii \sigma t)`, where :math:`\sigma` is
the angular frequency, this seperability can be exploited through the
following ansatzes:

#. The Eulerian perturbations to any scalar quantity :math:`f` and
   vector quantity :math:`\va` take the forms

   .. math::
      :label: e:osc-sol-form

      \begin{aligned}
      f'(r,\theta,\phi;t) &= \tf'(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t), \\
      \va'(r,\theta,\phi;t) &= \left[ \tar'(r) \, \ver + \tah'(r) \, r \, \nablah \right] Y^{m}_{\ell}(\theta,\phi) \, \exp(\ii \sigma t),
      \end{aligned}

   respectively. Here, the functions :math:`\tf'`, :math:`\tar'` and
   :math:`\tah'` describe the radial dependence of the perturbations;
   :math:`Y^{m}_{\ell}` is a spherical harmonic with harmonic degree
   :math:`\ell` and azimuthal order :math:`m`;

   .. math::

      \nablah \equiv \frac{1}{r} \left[ \vet \pderiv{}{\theta} + \frac{\vep}{\sin\theta} \pderiv{}{\phi} \right]

   is the horizontal part of the gradient operator; and
   :math:`(\ver,\vet,\vep)` are the unit basis vectors in the
   :math:`(r,\theta,\phi)` directions. 

#. The velocity perturbation vector can be expressed as

   .. math::

      \vv' = \pderiv{\vxi}{t},

   where the displacement perturbation vector :math:`\vxi` describes
   the spatial displacement of a fluid element from its equilibrium
   position. As a specific application of
   eqn. (:eq:`e:osc-sol-form`), :math:`\vxi` takes the form

   .. math::
      :label: e:osc-sol-form-vxi

      \vxi(r,\theta,\phi;t) = \left[ \txir(r) \, \ver + \txih(r) \, r \, \nablah \right] Y^{m}_{\ell}(\theta,\phi) \, \exp(\ii \sigma t),

   where :math:`\txir` and :math:`\txih` describe the radial
   dependence of the radial and horizontal displacement perturbations.

#. Lagrangian perturbations to scalar and vector quantities can be
   derived by applying eqn. (:eq:`e:osc-eul-lag`) to the expressions
   in eqn. (:eq:`e:osc-sol-form`). Specifically, the Lagrangian
   perturbations to any scalar quantity :math:`f` and vector quantity
   :math:`\va` take the forms

   .. math::

      \begin{aligned}
      \delta f(r,\theta,\phi;t) &= \delta \tf(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t), \\
      \delta \va(r,\theta,\phi;t) &= \left[ \delta \tar(r) \, \ver + \delta \tah(r) \, r \, \nablah \right] \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t),
      \end{aligned}

   respectively. Assuming a spherically symmetric equilibrium state,
   the functions :math:`\delta \tf`, :math:`\delta \tar` and
   :math:`\delta \tah` are related to their Lagrangian counterparts
   by

   .. math::

      \begin{aligned}
      \delta \tf &= \tf' + \txir \deriv{f}{r}, \\
      \delta \tar &= \tar' + \txir \deriv{a_{r}}{r}, \\
      \delta \tah &= \tah' + \frac{\txih}{r} a_{r}.
      \end{aligned}

Substituting the solution forms dictated by these anzatzes into the
linearized equations leads to a system of ordinary differential
equations for the radial functions :math:`\txir`, :math:`\txih`,
:math:`\tP'` etc. The mechanical (mass and momentum conservation)
equations become

.. math::
   :label: e:osc-sep-cont

   \trho' + \frac{1}{r^{2}} \deriv{}{r} \left( \rho r^{2} \txir \right) - \frac{\lambda}{r} \rho \txih = 0,

.. math::
   :label: e:osc-sep-mom

   \begin{aligned}
   -\sigma^{2} \rho \txir &= - \deriv{\tP'}{r} - \trho' \deriv{\Phi}{r} - \rho \deriv{\tPhi'}{r}, \\
   -\sigma^{2} \rho \txih &= - \frac{1}{r} \left( \tP' + \rho \tPhi' \right).
   \end{aligned}

Likewise, Poisson's equation becomes

.. math::
   :label: e:osc-sep-poisson

   \frac{1}{r^{2}} \deriv{}{r} \left( r^{2} \deriv{\tPhi'}{r} \right) - \frac{\lambda}{r^{2}} \tPhi' = 4 \pi G \trho'

and the heat equation becomes

.. math::

   -\ii \sigma T \delta \tS =
   \delta \tepsnuc
   - \deriv{\delta \tLrad}{M_{r}} +
   \lambda \frac{\txih}{r} \deriv{\Lrad}{M_{r}} +
   \frac{\lambda}{\rho r} \left( \delta \tFradh - \frac{\txih}{r} \right),

where

.. math::

   \delta \tLrad \equiv 4 \pi r^{2} \left( \delta \tFradr + 2 \frac{\txir}{r} \Fradr \right)

is the Lagrangian perturbation to the radiative luminosity. The radiative diffusion equation becomes

.. math::

   \begin{aligned}
   \delta\tFradr &= \left[
   4 \frac{\delta \tT}{T} - \frac{\delta\tkappa}{\kappa} + 2 \frac{\txir}{r} - \lambda \frac{\txih}{r} +
   \frac{\sderiv{(\delta \tT/T)}{\ln r}}{\sderiv{\ln T}{\ln r}} \right] \Fradr, \\
   \delta\tFradh &= \left[ \frac{1}{\sderiv{\ln T}{\ln r}} \frac{\delta \tT}{T} - \frac{\txir}{r} + \frac{\txih}{r} \right] \Fradr.
   \end{aligned}

Finally, the thermodynamic, nuclear and opacity relations become

.. math::

   \frac{\delta \trho}{\rho} = \frac{1}{\Gammi} \frac{\delta \tP}{P} - \upsT \frac{\delta \tS}{\cP},
   \qquad
   \frac{\delta \tT}{T} = \nabla_{\rm ad} \frac{\delta \tP}{P} + \frac{\delta \tS}{\cP},

.. math::

   \frac{\delta \tepsnuc}{\epsnuc} = \epsnucad \frac{\delta \tP}{P} + \epsnucS \frac{\delta \tS}{\cP},
   \qquad
   \frac{\delta \tkappa}{\kappa} = \kapad \frac{\delta \tP}{P} + \kapS \frac{\delta \tS}{\cP}.

In these equations,

.. math::
   :label: e:lambda-norot

   \lambda = \ell(\ell+1)

is the eigenvalue of the angular parts of the oscillation equations,
which enters here into the radial parts as a separation constant. It
is related to the local horizontal wavenumber by

.. math::

   k_{\rm h}^{2} = \frac{\lambda}{r^{2}}.
