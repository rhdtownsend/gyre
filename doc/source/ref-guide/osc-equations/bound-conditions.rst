.. _osc-bound-conds:

Boundary Conditions
===================

To form a closed system, the :ref:`separated equations <osc-sep-eqns>`
are augmented by algebraic relations at the inner and outer boundaries
of the computational domain, and possibly at interior points as well.

Inner Boundary
--------------

When the inner boundary is placed at the stellar origin (:math:`r=0`),
the requirement that solutions remain finite leads to the set of
regularity conditions

.. math::

   \begin{align}
   \txir - \ell \txih = 0, \\
   \ell \tPhi' - r \deriv{\tPhi'}{r} = 0, \\
   \delta \tS = 0.
   \end{align}
   
Sometimes it's desirable that the inner boundary is instead placed at
:math:`r > 0` --- for instance, to excise the stellar core from
the oscillation calculations. Then, there is much more flexibility in the
choice of inner boundary condition. Possible options include setting
:math:`\txir = 0` or :math:`\txih=0` instead of the first equation
above.

Outer Boundary
--------------

The outer boundary typically corresponds to the stellar surface. Under
the assumption that the density vanishes on and above this surface, the
gravitational potential must match onto a solution to Laplace's
equation that remains finite at infinity, leading to the potential
boundary condition

.. math::

   (\ell + 1) \tPhi' + r \deriv{\tPhi'}{r} = 0.

Likewise, the assumption that there is no external pressure acting on
the star (consistent with the vanishing surface density) gives the
momentum boundary condition

.. math::

   \delta \tP = 0.

Finally, the thermal boundary condition can be derived from the
equation 

.. math::

   T^{4}(\tau) = \frac{4\sigma}{3} \Fradr \left( \tau + \frac{2}{3} \right)

describing the thermal structure of an atmosphere under the combined
plane-parallel, grey, Eddington, local thermodynamic equilibrium and
radiative equilibrium approximations. Here, :math:`\tau` is the
vertical optical depth and :math:`\sigsb` the Stefan-Boltzmann
constant. Setting :math:`\tau=0` (again, consistent with the vanishing
surface density) and perturbing this equation yields the desired
boundary condition

.. math::

   4 \frac{\delta \tT}{T} = \frac{\delta \tFradr}{\Fradr}.

Complications arise when realistic stellar models are considered,
because these typically extend only out to an optical depth
:math:`\tau=2/3` (or some similar value) corresponding to the
photosphere. In such cases the density does not vanish at the nominal
stellar surface, and the outer boundary conditions must be modified to
account for the effects of the atmosphere region lying above the
surface. Many stellar oscillation codes, including GYRE, can adopt
more sophisticated formulations for the momentum boundary condition
--- for instance, based on the assumption that the outer atmosphere
has an isothermal stratification. However, the atmospheric effects on
the potential and thermal boundary conditions are usually neglected.

Internal Boundaries
-------------------

Internal boundaries arise when the density and other related
quantities show a radial discontinuity within the star. Across such a
discontinuity :math:`\txir`, :math:`\delta \tP` and :math:`\delta
\tFradr` must remain continuous\ [#continuous]_. Internal boundary
conditions on other perturbations follow from integrating the
:ref:`separated equations <osc-sep-eqns>` across the discontinuity,
resulting in

.. math::

   \begin{align}
   \tP^{\prime +} - \tP^{\prime -} &= \deriv{\Phi}{r} \left( \rho^{+} - \rho^{-} \right) \txir, \\
   \left. \deriv{\tPhi'}{r} \right|^{+} - \left. \deriv{\tPhi'}{r} \right|^{-} &= - 4 \pi G \left( \rho^{+} - \rho^{-} \right) \txir, \\
   \tT^{\prime +} - \tT^{\prime -} &= 0.
   \end{align}
   
Here, + (-) superscripts indicate quantities evaluated on the inner
(outer) side of the discontinuity.

.. rubric:: Footnotes

.. [#continuous] This is to ensure that the fluid doesn't 'tear', and
                 that pressure forces and radiative heating remain
                 finite.
