.. _osc-eqs:

*********************
Oscillation Equations
*********************

This chapter dervies the oscillation equations that GYRE solves.

.. toctree::
   :maxdepth: 2

.. _osc-eqs-fluid:

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
velocity; :math:`\Phi` is the gravitational potential; :math:`\epsnuc`
is the specific nuclear energy generation rate; and :math:`\vFrad` and
:math:`\vFcon` are the radiative and convectiveq energy fluxes. An
explicit expression for the radiative flux is provided by the
radiative diffusion equation,

.. math::

   \vFrad = \frac{c}{3\kappa\rho} \nabla (a T^{4}),

where :math:`\kappa` is the opacity and :math:`a` the radiation
constant.

Thermodynamic Relations
=======================

The fluid equations are augmented by the thermodynamic relationships
between the four state variables (:math:`P`, :math:`T`, :math:`\rho`
and :math:`S`). Only two of these are required to uniquely specify the
state (we assume that the composition remains fixed over an
oscillation cycle). In GYRE, :math:`p` and :math:`S` are adopted as
these primary variables, and the other two are presumed to be
derivable from them:

.. math::

   \rho = \rho(P, S), \qquad
   T = T(P, S).

Equilibrium State
=================

In a static equilibrium state the fluid velocity vanishes. The
momentum equation then becomes the hydrostatic equilibrium equation

.. math::

  \nabla P = - \rho \nabla \Phi.

Assume the equilibrium state is spherically symmetric, this simplifies to

.. math::

  \deriv{P}{r} = - \rho \deriv{\Phi}{r}.

Poisson's equation can be integrated once to yield

.. math::

   \deriv{\Phi}{r} = \frac{G}{r^{2}} \int 4 \pi \rho r^{2} \, \diff{r} = \frac{G M_{r}}{r^{2}},

and so the hydrostatic equilibrium equation becomes

.. math::

  \deriv{P}{r} = - \rho \frac{G M_{r}}{r^{2}}.

The heat equation in the equilibrium state is

.. math::

   \rho T \pderiv{S}{t} = \rho \epsnuc - \nabla \cdot (\vFrad + \vFcon).

If the star is in thermal equilibrium then the left-hand side
vanishes, and the nuclear heating rate balances the flux divergence
term.

.. _osc-eqs-linearize:

Linearized Equations
====================

Applying an Eulerian (fixed position, denoted by a prime) perturbation
to the mass and momentum conservation equations, they linearize about
the static equilibrium state as

.. math::

   \rho' + \nabla \cdot ( \rho \vv' ) = 0,

.. math::

   \rho \pderiv{\vv'}{t} = - \nabla P' + \frac{\rho'}{\rho} \nabla P - \rho \nabla \Phi'.

(in these expressions, the absence of a prime now denotes an
equilibrium quantity).  Likewise, Poisson's equation becomes

.. math::

   \nabla^{2} \Phi' = 4 \pi G \rho'

Applying a Lagrangian (fixed mass element, denoted by a
:math:`\delta`) perturbation to the heat equation, and neglecting[#freeze]_ the
convective heating term :math:`\delta (\rho^{-1} \nabla \cdot
\vFcon)`, it linearizes about the equilibrium state as

.. math::

   T \pderiv{\delta S}{t} = \epsad \frac{\delta P}{P} + \epsS \frac{\delta S}{c_{P}} - \delta \left( \frac{1}{\rho} \nabla \cdot \vFrad \right).

The radiative diffusion equation likewise linearizes as

.. math::

in the radial direction, and

.. math::

in the horizontal direction (here, t

The thermodynamic relations likewise linearize to 

.. math::

   \frac{\delta \rho}{\rho} = \frac{1}{\Gamma_{1}} \frac{\delta P}{P} - \upsilon_{T} \frac{\delta S}{c_{P}},

.. math::

   \frac{\delta T}{T} = \nabla_{\rm ad} \frac{\delta P}{P} + \frac{\delta S}{c_{P}}.

The thermodynamic partial derivatives are defined as

.. math::

   \Gamma_{1} = \left( \pderiv{\ln P}{\ln \rho} \right)_{S} \qquad
   \upsilon_{T} = \left( \pderiv{\ln \rho}{\ln T} \right)_{P} \qquad
   c_{P} = \left( \pderiv{S}{\ln T} \right)_{P} \qquad
   \nabla_{\rm ad} = \left( \pderiv{\ln T}{\ln P} \right)_{S},

and the nuclear and opacity partials are

.. math::

   \epsad = \left( \pderiv{\ln \epsnuc}{\ln P} \right)_{\rm ad} \qquad
   \epsS = \left( \pderiv{\ln \epsnuc}{S} \right)_{P} \qquad
   \kapad = \left( \pderiv{\ln \kappa}{\ln P} \right)_{\rm ad} \qquad
   \kapS = \left( \pderiv{\ln \kappa}{S} \right)_{P} \qquad
   
.. _osc-eqs-sep:

Separation
==========

With a separation of variables in spherical-polar coordinates
:math:`(r,\theta,\phi)`, and assuming an oscillatory time (:math:`t`)
dependence with angular frequency :math:`\sigma`, solutions to the
linearized fluid equation can be expressed as

.. math::

   \xir(r,\theta,\phi;t) = \operatorname{Re} \left[ \sqrt{4\pi} \, \txir(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t) \right],

.. math::

   \vxih(r,\theta,\phi;t) = \operatorname{Re} \left[ \sqrt{4\pi} \, \txih(r) \, r \nablah Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t) \right],

.. math::

   f'(r,\theta,\phi;t) = \operatorname{Re} \left[ \sqrt{4\pi} \, \tf'(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t) \right]

Here, :math:`\xir` is the radial component of the displacement
perturbation vector :math:`\vxi`, and :math:`\vxih` is the
corresponding horizontal (polar and azimuthal) part of this vector;
:math:`\nablah` is the horizontal part of the spherical-polar gradient
operator; :math:`Y^{m}_{\ell}` is the spherical harmonic with harmonic
degree :math:`\ell` and azimuthal order :math:`m`; and :math:`f`
stands for any perturbable scalar. The displacement perturbation
vector is related to the velocity perturbation via

.. math::

   \vv' = \pderiv{\vxi}{t}

Substituting the above solution forms into the linearized
equations, the mechanical (mass and momentum conservation) equations
become

.. math::

   \trho' + \frac{1}{r^{2}} \deriv{}{r} \left( r^{2} \txir \right) - \frac{\ell(\ell+1)}{r} \rho \txih = 0,

.. math::

   -\sigma^{2} \rho \txir = - \deriv{\tP'}{r} + \frac{\trho'}{\rho} \deriv{P}{r} - \rho \deriv{\tPhi'}{r},

.. math::

   -\sigma^{2} \rho r \txih = - \tP' - \rho \tPhi'.

Likewise, Poisson's equation becomes

.. math::

   \frac{1}{r^{2}} \deriv{}{r} \left( r^{2} \deriv{\tPhi'}{r} \right) - \frac{\ell(\ell+1)}{r^{2}} \txih = 4 \pi G \trho'

and the thermal (heat and thermodynamic) equations become

.. math::

   -\ii \sigma T \Delta \tS = \epsad \frac{\delta \tP}{P} + \epsS \frac{\delta \tS}{c_{P}} - \delta \left( \frac{1}{\rho} \nabla \cdot \vF \right),

      T \pderiv{\delta S}{t} = \delta \epsnuc - \delta \left( \frac{1}{\rho} \nabla \cdot \vF \right),

.. math::

   \frac{\delta \rho}{\rho} = \frac{1}{\Gamma_{1}} \frac{\delta P}{P} - \upsilon_{T} \frac{\delta S}{c_{P}},

.. math::

   \frac{\delta T}{T} = \nabla_{\rm ad} \frac{\delta P}{P} + \frac{\delta S}{c_{P}}.

.. _osc-eqs-dimless:

Dimensionless Form
==================

TBD

.. rubric:: Footnotes

.. [#freeze] This is known as the *frozen convection*
             approximation. GYRE offers multiple ways to freeze
             convection; the one here is the default.
