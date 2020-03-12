.. _oscillation-equations:

*****************************
Mathmatical Formalism of GYRE
*****************************

This chapter details the mathematical formalism oscillation on which
GYRE is built.

.. toc::
   :maxdepth: 2

Fluid Equations
===============

The basis for most subsequent equations are the fluid conservation
equations, comprising the continuity equation

.. math::

   \left( \pderiv{}{t}{} + \vv \cdot \nabla \right) \rho = - \rho \nabla \cdot \vv,

the momentum equation

.. math::

   \rho \left( \pderiv{}{t}{} + \vv \cdot \nabla \right) \vv = - \nabla p - \rho \nabla \Phi,

and the heat equation 

.. math::

   \rho T \left( \pderiv{}{t} + \vv \cdot \nabla \right) S = \rho \epsnuc - \nabla \cdot \vF.

Here, :math:`\rho`, :math:`p`, :math:`T`, :math:`S` and :math:`\vv`
are the fluid density, pressure, temperature, specific entropy and
velocity; while :math:`\Phi` is the gravitational potential,
:math:`\epsnuc` is the specific nuclear energy generation rate and
:math:`\vF` the energy flux. The gravitational potential itself
obeys Poisson's equation

.. math::

   \nabla^{2} \Phi = 4 \pi G \rho.

The energy flux is the sum of the radiative (:math:`\vFrad`) and convective (:math:`\vFcon`) fluxes,

.. math::

  \vF = \vFrad + \vFcon;

the radiative flux is given by the radiative diffusion equation,

.. math::

   \vFrad = \frac{c}{3\kappa\rho} \nabla (a T^{4})
   
Equilibrium State
=================

In a static equilibrium state the fluid velocity vanishes. The
momentum equation then becomes the hydrostatic equilibrium equation

.. math::

   \nabla P + \rho \nabla \Phi = 0.

Linearization
=============

Linearizing the fluid equations about an equilibrium state




