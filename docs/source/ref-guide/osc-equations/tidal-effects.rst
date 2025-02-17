.. _osc-tidal:

Tidal Effects
=============

To simulate the effects of tidal forcing by a companion, the
:program:`gyre_tides` frontend solves a modified form of the
linearized momentum equation (:eq:`e:osc-lin-mom`), namely

.. math::

   \rho \pderiv{\vv'}{t} = - \nabla P' - \rho' \nabla P - \rho \nabla \Phi' - \rho \nabla \PhiT.

The final term on the right-hand side represents the external force
arising from the tidal gravitational potential :math:`\PhiT`.

Tidal Potential
---------------

The tidal potential can be expressed via the superposition

.. math::
   :label: e:tidal-pot

   \PhiT = \sum_{\ell=2}^{\infty} \sum_{m=-\ell}^{\ell} \sum_{k=-\infty}^{\infty} \PhiTlmk.

of partial tidal potentials defined by

.. math::

  \PhiTlmk \equiv
  - \epsT \,
  \frac{G\Mstar}{\Rstar} \,
  \cbar_{\ell,m,k}
  \left( \frac{r}{\Rstar} \right)^{\ell} Y^{m}_{\ell}(\theta, \phi) \,
  \exp(- \ii k \Oorb t).

(the summation over :math:`\ell` and :math:`m` comes from a multipolar
space expansion of the potential, and the summation over :math:`k`
from a Fourier time expansion). Here,
   
.. math::

   \epsT = \left( \frac{\Rstar}{a} \right)^{3} q = \frac{\Oorb^{2} \Rstar^{3}}{G\Mstar} \frac{q}{1+q}

quantifies the overall strength of the tidal forcing, in terms of the
companion's mass :math:`q M`, semi-major axis :math:`a` and orbital
angular frequency :math:`\Oorb`. These expressions, and the definition
of the tidal expansion coefficients :math:`\cbar_{\ell,m,k}`, are presented in
greater detail in :ads_citet:`sun:2023`.

Separated Equations
-------------------

Because the tidal potential (:eq:`e:tidal-pot`) superposes many
different spherical harmonics, the solution forms
(:eq:`e:osc-sol-forms`) must be replaced by the more-general
expressions

.. math::
   :label: e:tidal-sol-forms

   \begin{aligned}
   \xir(r,\theta,\phi;t) &= \sum_{\ell,m,k} \txirlmk(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii k \Oorb t), \\
   \xit(r,\theta,\phi;t) &= \sum_{\ell,m,k} \txihlmk(r) \, \pderiv{}{\theta} Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii k \Oorb t), \\
   \xip(r,\theta,\phi;t) &= \sum_{\ell,m,k} \txihlmk(r) \, \frac{\ii m}{\sin\theta} Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii k \Oorb t), \\
   f'(r,\theta,\phi;t) &= \sum_{\ell,m,k} \tflmk'(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii k \Oorb t)
   \end{aligned}

(the notation for the sums has been abbreviated). Substituting these
solution forms into the :ref:`linearized equations <osc-linear-eqns>`,
and taking advantage of the orthonormality of the spherical harmonics,
leads to a fully separated set of differential equations for each
combination of :math:`\ell`, :math:`m` and :math:`k`. A given set
resembles the regular :ref:`separated equations <osc-sep-eqns>`, with
just a couple changes:

- The perturbation :math:`\tPhi'` is replaced by :math:`\tPsi' \equiv
  \tPhi' + \tPhiT`, representing the total (self + tidal)
  gravitational potential perturbation.
- Rather than being an eigenvalue parameter, the oscillation frequency
  is set by :math:`\sigma = k \Oorb`, representing the forcing
  frequency of the partial tidal potential in an inertial frame.

The latter change means that the dimensionless frequency (:eq:`e:omega`) becomes

.. math::

   \omega = \alphafrq \, k \Oorb \sqrt{\frac{\Rstar^{3}}{G\Mstar}},

where :math:`\alphafrq` is an additional term introduced to allow
tuning of the tidal forcing frequency (see the :nml_n:`alpha_frq` parameter
in the :ref:`tidal-params` section).

Boundary Conditions
-------------------

The boundary conditions accompanying the separated equations for a
given :math:`\{\ell,m,k\}` combination resemble those :ref:`presented
previously <osc-bound-conds>`, except that the outer potential
boundary condition becomes

.. math::

   (\ell + 1) \tPsi' + r \deriv{\tPsi'}{r} = (2 \ell + 1) \tPhiTlmk,

where

.. math::
   :label: e:tidal-part-pot
   
   \tPhiTlmk \equiv - \epsT \,
   \frac{G\Mstar}{\Rstar} \,
   \cbar_{\ell,m,k}
   \left( \frac{r}{\Rstar} \right)^{\ell}.

describes the radial dependence of the partial tidal potential.
