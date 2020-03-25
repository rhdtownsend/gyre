.. _analytic-solution:

Analytic Solution
=================

The :ref:`stretched-string BVP <stretched-string>` is straithgforward
to solve analytically. General solutions of the wave equation take the
form of traveling waves,

.. math::

  y(x,t) = A \exp [\ii (k x - \sigma t) ],

where :math:`A` an arbitrary constant, and the frequency
:math:`\sigma` and wavenumber :math:`k` are linked by the dispersion
relation

.. math::

  \sigma^{2} = c^{2} k^{2}.

The phase velocity of these waves is :math:`\sigma/k = \pm c`.

To satisfy the boundary condition at :math:`x=0`, we combine
traveling-wave solutions with opposite-sign wavenumbers

.. math::

   y(x,t) = A \exp [\ii (k x - \sigma t) ] - A \exp [\ii (- k x - \sigma t) ] = B \sin(k x) \exp ( - \ii \sigma t),

where :math:`B = 2A`. For the boundary condition at :math:`x=L` to be
satisfied simultaneously,

.. math::

   \sin(k L) = 0,

and so

.. math::

   k L = n \pi

where :math:`n` is a non-zero integer (we exclude :math:`n=0` because
it leads to the trivial solution :math:`y(x,t)=0`). . Combining this with
the dispersion relation, we find the normal-mode eigenfrequencies of
the stretched-string BVP are

.. math::
   :label: analytic-eigenfreqs

   \sigma = n \frac{\pi c}{L},

and the corresponding eigenfunctions are

.. math::
   :label: analytic-eigenfuncs

   y(x,t) = B \sin \left( \frac{n \pi x}{L} \right) \exp ( - \ii \sigma t).

The index :math:`n` uniquely labels the modes, and :math:`y(x,t)`
exhibits :math:`n-1` nodes in the open interval :math:`x=(0,L)`.
