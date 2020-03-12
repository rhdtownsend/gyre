.. _stretched-string:

The Stretched String Problem
============================

We'll start our discussion of GYRE by considering the analogous (but
much simpler) problem of finding normal-mode eigenfrequencies and
eigenfunctions for waves on a stretched string clamped at both
ends. Let the string have mass per unit length :math:`\rho` and
tension :math:`T`; then, the wave equation describing the transverse
string displacement :math:`y(x,t)` at spatial position :math:`x` and
time :math:`t` is

.. math::

   \npderiv{y}{x}{2} = - \frac{1}{c^{2}} \npderiv{y}{t}{2},

with :math:`c \equiv (T/\rho)^{1/2}`. If the string is clamped at
:math:`x=0` and :math:`x=L`, then the wave equation together with the boundary conditions

.. math::
   y(0,t) = 0 \qquad
   y(L,t) = 0

comprise a two-point boundary value problem (BVP).
