.. _numerical-solution:

Numerical Solution
==================

Now let's see how we might go about solving the :ref:`stretched-string
BVP <stretched-string>` numerically, using an approach very similar to
the one implemented in GYRE's for solving the oscillation equations.

Separation
----------

We begin by performing a separation of variables on the wave equation,
assuming trial solutions of the form

.. math::

   y(x,t) = \tilde{y}(x) \exp (-\ii \sigma t),

where :math:`\tilde{y}(x)` is some function of :math:`x` alone. Then,
the wave equation reduces to an ordinary differential equation (ODE)
for :math:`\tilde{y}`,

.. math::

   \nderiv{\tilde{y}}{x}{2} = \frac{\sigma^{2}}{c^{2}} \tilde{y}.

.. _discretization:

Discretization
--------------

To solve the ODE, we discretize it to form a set of difference
equations. The discretization involves transforming the continuous
function :math:`\tilde{y}(x)` into a finite set of :math:`N` values
:math:`\{\tilde{y}_{1},\tilde{y}_{2},\ldots,\tilde{y}_{N}\}`, representing the function
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
-------------

To find solutions to the linear system, we first write it in matrix form as

.. math:: 
   :label: linear-sys

   \mS \vu = \mathbf{0},

where :math:`\vu` is the vector with components

.. math::

   \vu = 
   \begin{pmatrix}
   \tilde{y}_{1} \\
   \tilde{y}_{2} \\
   \vdots \\
   \tilde{y}_{N-1} \\
   \tilde{y}_{N}
  \end{pmatrix}

and the 'system matrix' :math:`\mS` is an :math:`N \times N` tridiagonal matrix
with components

.. math::

   \mS = 
   \begin{pmatrix}
   1 & 0 & 0 & \cdots & 0 & 0 & 0 \\
   1 & \sigma^{2} \tau^{2} - 2 & 1 & \cdots & 0 & 0 & 0 \\
   \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
   0 & 0 & 0 & \cdots & 1 & \sigma^{2} \tau^{2} - 2 & 1 \\
   0 & 0 & 0 & \cdots & 0 & 0 & 1
   \end{pmatrix}.

Here we've introduced 

.. math::

   \tau \equiv \frac{L}{N-1} c,

representing the sound crossing time of a single cell.

Equation :eq:`linear-sys` is a :wiki:`homogeneous linear system
<System_of_linear_equations#Homogeneous_systems>`, meaning that it
only has non-trivial solutions :math:`\vu` when the determinant of
:math:`\mS` vanishes. With this in mind, we formulate the
characteristic equation for the BVP,

.. math::

   \Dfunc(\sigma) = 0

where :math:`\Dfunc(\sigma) \equiv \det(\mS)` is a
discriminant function whose roots are the characteristic frequencies
(*eigenfrequencies*) of the stretched-string BVP.

.. _fig-discrim-func:

.. figure:: fig_discrim_func.svg
   :alt: Plot showing the discriminant function versus frequency
   :align: center

   Plot of the discriminant function :math:`\Dfunc(\sigma)` as a
   function of the frequency :math:`\sigma`, for the stretched-string BVP
   with :math:`N=50`. The orange dots highlight where
   :math:`\Dfunc=0`. The function has been scaled so that
   :math:`\Dfunc(0) = 1`. (:download:`Source
   <fig_discrim_func.py>`)

:numref:`fig-discrim-func` plots the discriminant function for the BVP
discretized on a spatial grid of :math:`N=50` points. The roots
(zeros) of the function are highlighted by the orange markers; they
fall very close to the values :math:`\sigma = \pi c/L, 2 \pi c/L,
\ldots` predicted by the :ref:`analytic-solution`.

Finding Eigenfrequencies
------------------------

While :numref:`fig-discrim-func` is useful for visalizing
:math:`\Dfunc`, it's not the best way to find
eigenfrequencies. Instead, we can rely on well-established techniques
for isolating and refining roots of monovariate functions.

First, we evaluate a finite set of :math:`M` values
:math:`\{\Dfunc_{1},\Dfunc_{2},\ldots,\Dfunc_{M}\}`, representing the
discriminant function sampled on the discrete frequency grid
:math:`\{\sigma_{1},\sigma_{2},\ldots,\sigma_{M}\}`. Then, we inspect
the signs of adjacent values :math:`(\Dfunc_{j},\Dfunc_{j+1})`. If
these differ, then we know that a root of the discriminant function
must lie in the interval :math:`(\sigma_{j},\sigma_{j+1})` --- we have
*bracketed* a root. :numref:`fig-discrim-brackets` demonstrates the
process of root bracketing for a frequency grid covering the plotted
frequency interval with :math:`M=32` uniformly spaced points; it
highlights five brackets containing the five roots shown previously in
:numref:`fig-discrim-func`.

.. _fig-discrim-brackets:

.. figure:: fig_discrim_brackets.svg
   :alt: Plot showing the discriminant function versus frequency, with root brackets indicated
   :align: center

   Plot of the discriminant values :math:`\{\Dfunc\}` on the discrete
   frequency grid :math:`\{\sigma\}`, for the stretched-string BVP
   with :math:`N=50` and :math:`M=32`. The orange halos indicate
   adjacent points that bracket a root
   :math:`\Dfunc=0`. (:download:`Source <fig_discrim_brackets.py>`)

Once a bracket is established for a given root, it can be narrowed
through a process of iterative refinement until the root is converged
upon. There are a variety of well-known root-finding algorithms that
perform this refinement; the :wiki:`bisection method` is conceptually
the simplest, but approaches such as :wiki:`Brent's method` can be
much more efficient. For the brackets plotted in
:numref:`fig-discrim-brackets`, :numref:`numerical-eigenfreqs` compares
the eigenfrequencies found using Python's
:py:func:`scipy.optimize.brentq` function, against the analytic values
predicted by equation :eq:`analytic-eigenfreqs`.

.. _numerical-eigenfreqs:

.. csv-table:: Numerical and analytic eigenfrequencies, in units of
   :math:`\pi c/L`, for the stretched-string BVP with
   :math:`N=50`. (:download:`Source <discrim_roots.py>`)
   :widths: 20 40 40
   :align: center
   :file: discrim_roots.csv

Eigenfunction Reconstruction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each of the eigenfrequencies found, we find the corresponding
eigenfunction by solving the linear system :eq:`linear-sys`. Because
:math:`\det(\mS)` is now zero, this system is guaranteed to have a
non-trivial solution. The solution vector :math:`\vu` resides in the
:wiki:`null space` of :math:`\mS`, and we can use standard numerical
techniques (e.g., :wiki:`singular value decomposition`) to evaluate
it. Then, the :math:`k`'th element of :math:`\vu` corresponds to the
eigenfunction sampled at the :math:`k`'th spatial grid point:

.. math::

   (\vu)_{k} = \tilde{y}_{k} \equiv \tilde{y}(x_{k})

.. _fig-eigenfuncs:

.. figure:: fig_eigenfuncs.svg
   :alt: Plot showing eigenfunctions for the first three modes
   :align: center

   Plot of the eigenfunctions :math:`\tilde{y}` as a function of
   spatial coordinate :math:`x`, for the first three modes of the
   stretched-string BVP with :math:`N=50`. The discrete points show
   the numerical functions, and the solid lines the corresponding
   analytic functions. (:download:`Source <fig_eigenfuncs.py>`)

:numref:`fig-eigenfuncs` plots the eigenfunctions found in this way
for the first three modes (:math:`n=1,\ldots,3`) given in
:numref:`numerical-eigenfreqs`. Also shown are the corresponding
analytic solutions given by equation :eq:`analytic-eigenfuncs`. The
agreement between the two is good.

