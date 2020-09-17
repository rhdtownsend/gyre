.. _equilib-state:

Equilibrium State
=================

In the static equilibrium state the fluid velocity :math:`\vv`
vanishes. The momentum equation then becomes the hydrostatic
equilibrium equation

.. math::

  \nabla P = - \rho \nabla \Phi.

Assuming the equilibrium is spherically symmetric, this simplifies to

.. math::

  \deriv{P}{r} = - \rho \deriv{\Phi}{r}.

Poisson's equation can be integrated once to yield

.. math::

   \deriv{\Phi}{r} = \frac{G}{r^{2}} \int 4 \pi \rho r^{2} \, \diff{r} = \frac{G M_{r}}{r^{2}},

where the second equality introduces the interior mass

.. math::

   M_{r} \equiv \int 4 \pi \rho r^{2} \, \diff{r}.

The hydrostatic equilibrium equation thus becomes

.. math::

  \deriv{P}{r} = - \rho \frac{G M_{r}}{r^{2}}.

The heat equation in the equilibrium state is

.. math::

   \rho T \pderiv{S}{t} = \rho \epsnuc - \nabla \cdot (\vFrad + \vFcon).

If the star is in thermal equilibrium then the left-hand side
vanishes, and the nuclear heating rate balances the flux divergence
term. Again assuming spherical symmetry, this is written

.. math::

   \deriv{}{r} \left( \Lrad + \Lcon \right) = 4 \pi r^{2} \rho \epsnuc,

where

.. math::

   \Lrad \equiv 4 \pi r^{2} \Fradr, \qquad
   \Lcon \equiv 4 \pi r^{2} \Fconr

are the radiative and convective luminosities, respectively.
