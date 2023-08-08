.. _osc-conv:

Convection Effects
==================

The oscillation equations presented in the preceding sections neglect
the thermal and mechanical effects of convection. GYRE provides
functionality for controlling how the thermal effects are suppressed,
and how the mechanical effects can be included in a limited way.

.. _osc-conv-frozen:

Frozen Convection
-----------------

In the derivation of the :ref:`linearized equations
<osc-linear-eqns>`, a term :math:`\delta (\rho^{-1} \nabla \cdot
\vFcon)` is dropped from the perturbed heat equation. This is known as
a *frozen convection* approximation, and is grounded in the assumption
that the energy transport by convection remains unaffected affected by
the pulsation. There's more than one way to freeze convection;
:ads_citet:`pesnell:1990` presents a systematic review of different
approaches. GYRE currently implements a subset of these:

* Pesnell's case 1, neglecting :math:`\delta (\rho^{-1} \nabla \cdot \vFcon)` in the perturbed heat equation.
* Pesnell's case 4, neglecting :math:`\delta \Lcon` (the Lagrangian
  perturbation to the convective luminosity) in the perturbed heat
  equation.

For further details, see the :nml_n:`conv_scheme` parameter in the
:ref:`osc-params` section.

.. _osc-conv-turb:

Turbulent Damping
-----------------

The Reynolds number in stars is very large, and thus convection tends
to be turbulent. Following the treatment by
:ads_citet:`savonije:2002`, GYRE can partially incorporate the
mechanical effects of this turbulence by adding a term

.. math::

   f_{r,{\rm visc}} = - \frac{1}{r^{2}} \pderiv{}{r} \left( \rho \nu r^{2} \pderiv{v'_{r}}{r} \right)

to the radial component of the linearized momentum equation
(:eq:`e:osc-lin-mom`), representing the viscous force arising from
radial fluid motion. Because this term depends on :math:`v'_{r}`, it
is phase-shifted by a quarter cycle relative to the other terms in the
equation, and acts like a drag force that damps oscillations. The
turbulent viscous coefficient :math:`\nu` is evaluated as

.. math::

   \nu = \frac{(\alphatrb H_{P})^{2}}{\tconv}
   \left[ 1 + \tconv \frac{\sigma}{2\pi} \right]^{-1},

where :math:`H_{P}` is the pressure scale height, :math:`\alphatrb` is
the turbulent mixing length (in units of :math:`H_{P}`), and
:math:`\tconv` the convective turnover timescale. This expression is
adapted from equation (18) of :ads_citet:`savonije:2002`, with an
exponent :math:`s=1`.

In GYRE :math:`\alphatrb` is implemented as a switch (see the
:ref:`osc-physics-switches` section). A reasonable choice is to set
this parameter equal to the MLT mixing length parameter
:math:`\alpha_{\rm MLT}` of the stellar model.
