.. _tidal-linear-eqns:

Linearized Equations
====================

The linearized tidal equations are similar to the :ref:`linearized
oscillation equations <osc-linear-eqns>`, but include an extra term in
the momentum equation :eq:`e:osc-lin-mom` representing the
tidal force exerted by the companion:

.. math::

   \rho \pderiv{\vv'}{t} = - \nabla P' - \rho' \nabla P - \rho \nabla \Phi' - \rho \nabla \PhiT.

The tidal potential :math:`\PhiT` is expressed as a superposition

.. math::
   :label: e:tidal-pot

   \PhiT = \sum_{\ell=2}^{\infty} \sum_{m=-\ell}^{\ell} \sum_{k=-\infty}^{\infty} \PhiTlmk,

of partial tidal potentials defined by

.. math::

  \PhiTlmk = 
  - \epsT \,
  \frac{GM}{R} \,
  \cbar_{\ell,m,k}
  \left( \frac{r}{R} \right)^{\ell} Y^{m}_{\ell}(\theta, \phi) \,
  \exp(- \ii k \Oorb t).

Here,
   
.. math::

   \epsT = \left( \frac{R}{a} \right)^{3} = \frac{\Oorb R^{3}}{GM} \frac{q}{1+q}

quantifies the overall strength of the tidal forcing, in terms of the
companion's mass :math:`q M`, semi-major axis :math:`a` and orbital
angular frequency :math:`\Oorb`. These expressions, and the definition
of the tidal expansion coefficients :math:`\cbar_{\ell,m,k}`, are presented in
greater detail in :ads_citet:`sun:2023`.
