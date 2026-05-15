.. _osc-conv:

.. nml:group:: osc
   :no-target:

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

For further details, see the :nml:option:`conv_scheme
<osc.conv_scheme>` option of the :nml:group:`osc` namelist group.

.. _osc-conv-turb:

Turbulent Damping
-----------------

The Reynolds number in stars is very large, and thus convection tends
to be turbulent. Following the treatment by
:ads_citet:`willems:2010`, GYRE can partially incorporate the
mechanical effects of this turbulence by adding a term

.. math::

   f_{r,{\rm visc}} = \frac{1}{r^{2}} \pderiv{}{r} \left( \rho \nu r^{2} \pderiv{v'_{r}}{r} \right)

to the radial component of the linearized momentum equation
(:eq:`e:osc-lin-mom`), representing the viscous force per unit volume
arising from radial fluid motions. Because this term depends on
:math:`v'_{r}`, it is phase-shifted by a quarter cycle relative to the
other terms in the equation, and acts like a drag force that damps
oscillations. The turbulent viscosity coefficient :math:`\nu` is
evaluated as

.. math::

   \nu = \frac{L^{2}}{\tconv} 
   \left[ 1 + \left( \tconv \frac{\sigma}{2\pi} \right)^{\alphacon} \right]^{-1},

where :math:`L` is the mixing length, and :math:`\tconv` is the local
convection turnover timescale. The term in square brackets acts to
reduce the viscosity when the tidal forcing occurs at a rate faster
than the turnover timescale. As discussed by
:ads_citet:`willems:2010`, different authors have proposed different
exponents :math:`\alphacon`; GYRE's default :math:`\alphacon=1` can be
over-ridden using the :nml:option:`alpha_con` option.

GYRE evaluates the mixing length as

.. math::

   L = \alphatrb \min(H_{P}, r),

where :math:`H_{P}` is the local pressure scale height, and
:math:`\alphatrb` is implemented as a switch (see the
:ref:`osc-physics-switches` section). A reasonable choice is to set
:math:`\alphatrb` equal to the MLT mixing length parameter
:math:`\alpha_{\rm MLT}` of the stellar model. To disable turbulent
damping completely, set :math:`\alphatrb` to zero (the default).

To estimate the convection turnover timescale, GYRE uses the simple formula

.. math::

   \tconv = \left[ \max\left(-N^{2}, 0\right) \right]^{-1/2}.
