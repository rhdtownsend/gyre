.. _linear-equations:

Linearized Equations
====================

Applying an Eulerian (fixed position, denoted by a prime) perturbation
to the mass and momentum conservation equations, they linearize about
the static equilibrium state as

.. math::

   \rho' + \nabla \cdot ( \rho \vv' ) = 0,

.. math::

   \rho \pderiv{\vv'}{t} = - \nabla P' - \rho' \nabla P - \rho \nabla \Phi'.

(in these expressions, the absence of a prime denotes an
equilibrium quantity).  Likewise, Poisson's equation becomes

.. math::

   \nabla^{2} \Phi' = 4 \pi G \rho'

Applying a Lagrangian (fixed mass element, denoted by a
:math:`\delta`) perturbation to the heat equation, and neglecting\ [#freeze]_ the
convective heating term :math:`\delta (\rho^{-1} \nabla \cdot
\vFcon)`, it linearizes about the equilibrium state as

.. math::

   T \pderiv{\delta S}{t} = \delta \epsnuc - 
   \frac{1}{\rho} \nabla \cdot \left[ \vFrad' + \vxi (\nabla \cdot \vFrad) \right].

Likewise applying an Eulerian perturbation to the radiative diffusion equation,

.. math::

   \vFrad' = \Fradr \left[
   \left(  - \frac{\kappa'}{\kappa} - \frac{\rho'}{\rho} + 4 \frac{T'}{T} \right) \ver
   + \frac{\nabla (T'/T)}{\sderiv{\ln T}{r}} \right]

where :math:`\ver` is the radial unit vector. The thermodynamic
relations linearize to

.. math::

   \frac{\delta \rho}{\rho} = \frac{1}{\Gammi} \frac{\delta P}{P} - \upsT \frac{\delta S}{\cP},
   \qquad
   \frac{\delta T}{T} = \nabad \frac{\delta P}{P} + \frac{\delta S}{\cP},

and the peturbations to the nuclear energy generation rate and opacity can be expressed as

.. math::

   \frac{\delta \epsnuc}{\epsnuc} = \epsad \frac{\delta P}{P} + \epsS \frac{\delta S}{\cP},
   \qquad
   \frac{\delta \kappa}{\kappa} = \kapad \frac{\delta P}{P} + \kapS \frac{\delta S}{\cP}.

In these expressions, Eulerian and Lagrangian perturbations to any
scalar quantity :math:`f` are related via

.. math::

   \frac{\delta f}{f} = \frac{f'}{f} + \frac{\xir}{r} \deriv{\ln f}{\ln r}.

Moreover, the thermodynamic partial derivatives are defined as

.. math::

   \Gammi = \left( \pderiv{\ln P}{\ln \rho} \right)_{S}, \quad
   \upsT = \left( \pderiv{\ln \rho}{\ln T} \right)_{P}, \quad
   \cP = \left( \pderiv{S}{\ln T} \right)_{P}, \quad
   \nabad = \left( \pderiv{\ln T}{\ln P} \right)_{S},

and the nuclear and opacity partials are

.. math::

   \epsad = \left( \pderiv{\ln \epsnuc}{\ln P} \right)_{\rm ad}, \quad
   \epsS = \cP \left( \pderiv{\ln \epsnuc}{S} \right)_{P}, \quad
   \kapad = \left( \pderiv{\ln \kappa}{\ln P} \right)_{\rm ad}, \quad
   \kapS = \cP \left( \pderiv{\ln \kappa}{S} \right)_{P}.

The latter can be calculated from corresponding density and
temperature partials via

.. math::

   \begin{gathered}
   \kapad = \frac{\kaprho}{\Gammi} + \nabad \kapT, \qquad
   \kapS = -\upsT \kaprho + \kapT, \\
   \epsad = \frac{\epsrho}{\Gammi} + \nabad \epsT, \qquad
   \epsS = -\upsT \epsrho + \epsT.
   \end{gathered}

.. rubric:: Footnotes

.. [#freeze] This is known as the *frozen convection*
             approximation. GYRE offers multiple ways to freeze
             convection (see the :ref:`osc-params` section); the one
             here is the default.
   
