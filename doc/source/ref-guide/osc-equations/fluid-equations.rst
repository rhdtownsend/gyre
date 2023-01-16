.. _osc-fluid-eqns:

Fluid Equations
===============

The starting point is the fluid equations, comprising the conservation
laws for mass

.. math::

   \pderiv{\rho}{t} + \cdot \nabla \left( \rho \vv \right) = 0

and momentum

.. math::

   \rho \left( \pderiv{}{t} + \vv \cdot \nabla \right) \vv = -\nabla P - \rho \nabla \Phi;

the heat equation

.. math::

   \rho T \left( \pderiv{}{t} + \vv \cdot \nabla \right) S = \rho \epsnuc - \nabla \cdot (\vFrad + \vFcon);

and Poisson's equation

.. math::

   \nabla^{2} \Phi = 4 \pi G \rho.

Here, :math:`\rho`, :math:`P`, :math:`T`, :math:`S` and :math:`\vv`
are the fluid density, pressure, temperature, specific entropy and
velocity; :math:`\Phi` is the self-gravitational potential;
:math:`\epsnuc` is the specific nuclear energy generation rate; and
:math:`\vFrad` and :math:`\vFcon` are the radiative and convective
energy fluxes. An explicit expression for the radiative flux is
provided by the radiative diffusion equation,

.. math::

   \vFrad = - \frac{c}{3\kappa\rho} \nabla (a T^{4}),

where :math:`\kappa` is the opacity and :math:`a` the radiation
constant.

The fluid equations are augmented by the thermodynamic relationships
between the four state variables (:math:`P`, :math:`T`, :math:`\rho`
and :math:`S`). Only two of these are required to uniquely specify the
state (we assume that the composition remains fixed over an
oscillation cycle). In GYRE, :math:`P` and :math:`S` are adopted as
these primary variables\ [#choice]_, and the other two are presumed to be
derivable from them:

.. math::

   \rho = \rho(P, S), \qquad
   T = T(P, S).

The nuclear energy generation rate and opacity are likewise presumed
to be functions of the pressure and entropy:

.. math::

   \epsnuc = \epsnuc(P, S), \qquad
   \kappa = \kappa(P, S).
   
.. rubric:: Footnotes

.. [#choice] This may seem like a strange choice, but it simplifies
             the switch between adiabatic and non-adiabatic
             calculations
