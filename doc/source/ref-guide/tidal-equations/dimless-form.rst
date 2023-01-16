.. _tidal-dimless-form:

Dimensionless Formulation
=========================

The dimensionless formulation of the tidal equations is almost
identical to the :ref:`corresponding formulation <osc-dimless-form>`
of the :ref:`oscillation equations <osc-sep-eqns>`, differing only
in the definition of a couple of variables and a single boundary
condition.

Variables
---------

The definitions of the :math:`y_{3}` and :math:`y_{4}` dependent variables,
given in eqn. (:eq:`e:dimless`), are replaced by

.. math::

   \begin{align}
   y_{3} &= x^{2-\ell}\, \frac{\tPsi'}{gr}, \\
   y_{4} &= x^{2-\ell}\, \frac{1}{g} \deriv{\tPsi'}{r}.
   \end{align}

Boundary Conditions
-------------------

The outer potential boundary condition, the second line of
eqn. (:eq:`e:outer-bc`), is replaced by

.. math::

   \alphagrv U y_{1} + (\alphagrv \ell + 1) y_{3} + \alphagrv y_{4} = (2\ell+1) \yT,

where

.. math::
   :label: e:y_T

   \yT \equiv x^{2 - \ell} \frac{\tPhiTlmk}{gr}.
