.. _gyre-fundamentals:

*****************
GYRE Fundamentals
*****************

In this Chapter, we'll take a look at the fundamentals of
GYRE. Although it aims to be user-friendly, GYRE is nevertheless a
complex piece of software; thus, getting it to produce the 'best'
results requires some degree of insight into the algorithms it uses to
caclulate mode eigenfrequencies and eigenfunctions.

The Stretched String Problem
============================

We'll start our discussion by considering the analogous (but much
simpler) problem of finding normal-mode eigenfrequencies and
eigenfunctions for waves on a stretched string clamped at both
ends. Let the string have mass per unit length :math:`\rho` and
tension :math:`T`; then, the wave equation describing the transverse
string displacement :math:`y(x,t)` at spatial position :math:`x` and
time :math:`t` is

.. math::

   \pderiv{y}{x}{2} = - \frac{1}{c^{2}} \pderiv{y}{t}{2},

with :math:`c \equiv (T/\rho)^{1/2}`. If the string is clamped at
:math:`x=0` and :math:`x=L`, then the wave equation together with the boundary conditions

.. math::
   y(0,t) = 0 \qquad
   y(L,t) = 0

comprise a two-point boundary value problem (BVP).

Analytic Solution
-----------------

The BVP is straithgforward to solve analytically. General solutions of
the wave equation take the form of traveling waves,

.. math::

  y(x,t) = A \exp [\ii (k x - \sigma t) ],

where :math:`A` an arbitrary constant, and the frequency
:math:`\sigma` and wavenumber :math:`k` are linked by the dispersion
relation

.. math::

  \omega^{2} = c^{2} k^{2}.

The phase velocity of these waves is :math:`\omega/k = \pm c`.

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

where :math:`n` is an integer. Combining this with the dispersion
relation, we find the normal-mode eigenfrequencies of the
stretched-string BVP are

.. math::

   \sigma = n \frac{\pi c}{L}.

Numerical Solution
------------------

Now let's see how we might go about solving the same BVP numerically,
using an approach very similar to GYRE.

Separation
~~~~~~~~~~

We perform a separation of variables on the wave equation by assuming
trial solutions of the form

.. math::

   y(x,t) = \tilde{y}(x) \exp (-\ii \sigma t),

where :math:`\tilde{y}(x)` is some function of :math:`x` alone. Then,
the wave equation reduces to an ordinary differential equation (ODE)
for :math:`\tilde{y}`,

.. math::

   \deriv{\tilde{y}}{x}{2} = \frac{\sigma^{2}}{c^{2}} \tilde{y}.

Discretization
^^^^^^^^^^^^^^

To solve the ODE, we discretize it to form a set of difference
equations. The discretization involves transforming the continuous
function :math:`\tilde{y}(x)` into a finite set of :math:`N` values
:math:`\{y_{1},y_{2},\ldots,y_{N}\}`, representing the function
sampled on the discrete spatial grid
:math:`\{x_{1},x_{2},\ldots,x_{N}\}`.

For simplicity let's assume the grid is uniform, so that

.. math::

   x_{k+1} - x_{k} = \Delta x \equiv \frac{L}{N-1}
   \qquad (1 \leq k \leq N-1)

Then, the second derivative of :math:`\tilde{y}` can be approximated (to second order in :math:`\Delta x`) as

.. math::

   \left. \deriv{\tilde{y}}{x}{2} \right|_{x=x_{k}} \approx \frac{\tilde{y}_{k+1} - 2 \tilde{y}_{k} + \tilde{y}_{k-1}}{\Delta x^{2}}
   \qquad (2 \leq k \leq N-1).
   
This allows us to replace the ODE with :math:`N-2` difference
equations

.. math::

   \frac{\tilde{y}_{k+1} - 2 \tilde{y}_{k} + \tilde{y}_{k-1}}{\Delta x^{2}} = \frac{\sigma^{2}}{c^{2}} \tilde{y}_{k}
   \qquad (2 \leq k \leq N-1)

Taken together with the 2 boundary conditions

.. math::

   \tilde{y}_{1} = 0 \qquad
   \tilde{y}_{N} = 0,

we have a linear system of :math:`N` algebraic equations and :math:`N` unknowns.
   
Linear System
^^^^^^^^^^^^^

To find solutions to the linear system, we first write it in matrix form as

.. math::

   \mathbf{M} \mathbf{u} = \mathbf{0},

where :math:`\mathbf{u}` is the vector with components

.. math::

   \mathbf{u} = 
   \begin{pmatrix}
   \tilde{y}_{1} \\
   \tilde{y}_{2} \\
   \vdots \\
   \tilde{y}_{N-1} \\
   \tilde{y}_{N}
  \end{pmatrix}

and :math:`\mathbf{M}` is an :math:`N \times N` tridiagonal matrix
with components

.. math::

   \mathbf{M} = 
   \begin{pmatrix}
   1 & 0 & 0 & \cdots & 0 & 0 & 0 \\
   1 & \delta^{2} \omega^{2} - 2 & 1 & \cdots & 0 & 0 & 0 \\
   \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
   0 & 0 & 0 & \cdots & 1 & \delta^{2} \omega^{2} - 2 & 1 \\
   0 & 0 & 0 & \cdots & 0 & 0 & 1
   \end{pmatrix}.

Here, we've introduced the dimensionless frequency :math:`\omega` and grid-spacing parameter :math:`\delta` as

.. math::

   \omega \equiv \frac{\sigma L}{\pi c}, \qquad
   \delta \equiv \frac{\pi \Delta}{L} = \frac{\pi}{N-1}.


Root Finding
^^^^^^^^^^^^

.. math::

   \disc(\sigma) = 0

   




   

  

