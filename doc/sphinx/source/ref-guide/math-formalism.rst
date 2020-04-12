.. _math-formalism:

*********************
Mathmatical Formalism
*********************

This chapter details the mathematical formalism on which
GYRE is built.

.. toctree::
   :maxdepth: 2

.. _math-fluid:

Physical Formulation
====================

Fluid Equations
---------------

The basis for most subsequent equations are the fluid equations,
comprising the conservation laws for mass

.. math::

   \pderiv{\rho}{t} + \cdot \nabla \left( \rho \vv \right) = 0

and momentum

.. math::

   \rho \left( \pderiv{}{t} + \vv \cdot \nabla \right) \vv = -\nabla P - \rho \nabla \Phi,

the heat equation

.. math::

   \rho T \left( \pderiv{}{t} + \vv \cdot \nabla \right) S = \rho \epsnuc - \nabla \cdot \vF,

and Poisson's equation

.. math::

   \nabla^{2} \Phi = 4 \pi G \rho.

Here, :math:`\rho`, :math:`p`, :math:`T`, :math:`S` and :math:`\vv`
are the fluid density, pressure, temperature, specific entropy and
velocity; while :math:`\Phi` is the gravitational potential,
:math:`\epsnuc` is the specific nuclear energy generation rate and
:math:`\vF` the energy flux.
      
The energy flux is the sum of the radiative (:math:`\vFrad`) and
convective (:math:`\vFcon`) fluxes,

.. math::

  \vF = \vFrad + \vFcon;

the radiative flux is given by the radiative diffusion equation,

.. math::

   \vFrad = \frac{c}{3\kappa\rho} \nabla (a T^{4}),

where :math:`\kappa` is the opacity and :math:`a` the radiation
constant.

Equilibrium State
-----------------

In a static equilibrium state the fluid velocity vanishes. The
momentum equation then becomes the hydrostatic equilibrium equation

.. math::

  \nabla P + \rho \nabla \Phi = 0.

   
.. _math-linearize:

Linearized Equations
--------------------

Linearizing the fluid equations about the static equlibrium state,
they become

.. math::

   \rho' + \nabla \cdot ( \rho \vv' ) = 0,

.. math::

   \rho \pderiv{\vv'}{t} = - \nabla P' + \frac{\rho'}{\rho} \nabla P - \rho \nabla \Phi'

.. math::

   \rho T \pderiv{\delta S}{t} = \rho' \epsnuc + \rho \epsnuc' - \left( \nabla \cdot \vF \right)'.

.. math::

   \nabla^{2} \Phi' = 4 \pi G \rho'

.. math::

   \vFrad' = \left( \frac{c}{3\kappa\rho} \right)' \nabla (a T^{4}) + XXXXX

Here, a prime (:math:`'`) indicates an Eulerian (fixed position)
perturbation, and a :math:`\delta` indicates a Lagrangian (fixed mass
element) perturbation; the absence of either denotes an equilibrium
quantity. No :math:`\vv'` terms appear because the equilibrium state is
static.

.. _math-osc:

Separation
----------

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

Oscillation Equations
=====================

Dimensioned Form
~~~~~~~~~~~~~~~~

The oscillation equations follow from substituting the above solution
forms into the linearized equations:

.. math::

   \trho' + \frac{1}{r^{2}} \deriv{}{r} \left( r^{2} \txir \right) - \frac{\ell(\ell+1)}{r} \rho \txih = 0,

.. math::

   -\sigma^{2] \rho \txir = - \deriv{\tP'}{r} + \frac{\trho'}{\rho} \deriv{P}{r} - \rho \deriv{\tPhi'}{r},

.. math::

   -\sigma^{2} \rho r \thxi = - \tP' - \rho \tPhi',

