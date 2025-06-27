.. _numerical-string:

The Stretched String Problem
============================

We'll start our discussion of numerical methods by considering the
problem of finding normal-mode eigenfrequencies and eigenfunctions for
waves on a stretched string clamped at both ends. Let the string have
mass per unit length :math:`\rho` and tension :math:`T`; then, the
wave equation describing the transverse string displacement
:math:`y(x,t)` at spatial position :math:`x` and time :math:`t` is

.. math::

   \npderiv{y}{x}{2} = \frac{1}{c^{2}} \npderiv{y}{t}{2},

with :math:`c \equiv (T/\rho)^{1/2}`. If the string is clamped at
:math:`x=0` and :math:`x=L`, then the wave equation together with the boundary conditions

.. math::

   y(0,t) = 0 \qquad
   y(L,t) = 0

comprises a two-point boundary value problem (BVP).

.. _string-analytic:

Analytic Solution
-----------------

The stretched-string BVP is straightforward to solve
analytically. General solutions of the wave equation take the form of
traveling waves,

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
it corresponds to the trivial solution :math:`y(x,t)=0`). Combining
this with the dispersion relation, we find that the normal-mode
eigenfrequencies of the stretched-string BVP are

.. math::
   :label: analytic-eigenfreqs

   \sigma = n \frac{\pi c}{L},

and the corresponding eigenfunctions are

.. math::
   :label: analytic-eigenfuncs

   y_{n}(x,t) = B \sin \left( \frac{n \pi x}{L} \right) \exp ( - \ii \sigma t).

The index :math:`n` uniquely labels the modes, and :math:`y_{n}(x,t)`
exhibits :math:`n-1` nodes in the open interval :math:`x \in (0,L)`.

Separation
----------

Now let's see how we might go about solving the stretched-string BVP
numerically. We begin by performing a separation of variables on the
wave equation, assuming trial solutions of the form

.. math::
   :label: var-separation

   y(x;t) = \tilde{y}(x) \, \exp (-\ii \sigma t),

where :math:`\tilde{y}(x)` is a function of :math:`x` alone. Then,
the wave equation reduces to an ordinary differential equation (ODE)
for :math:`\tilde{y}`,

.. math::

   \nderiv{\tilde{y}}{x}{2} = - \frac{\sigma^{2}}{c^{2}} \tilde{y}.

Discretization
--------------

To solve the ODE, we discretize it to establish a set of difference
equations. The discretization involves transforming the continuous
function :math:`\tilde{y}(x)` into a finite set of :math:`N` values
:math:`\{\tilde{y}_{1},\tilde{y}_{2},\ldots,\tilde{y}_{N}\}`,
representing the function sampled on the discrete spatial grid
:math:`\{x_{1},x_{2},\ldots,x_{N}\}`.

For simplicity let's assume the grid is uniform, so that

.. math::

   x_{j+1} - x_{j} \equiv \Delta x = \frac{L}{N-1}
   \qquad (1 \leq j \leq N-1).

Then, the second derivative of :math:`\tilde{y}` can be approximated to second order in :math:`\Delta x` as

.. math::

   \left. \nderiv{\tilde{y}}{x}{2} \right|_{x=x_{j}} \approx \frac{\tilde{y}_{j+1} - 2 \tilde{y}_{j} + \tilde{y}_{j-1}}{\Delta x^{2}}
   \qquad (2 \leq j \leq N-1).
   
This allows us to replace the ODE with :math:`N-2` difference
equations

.. math::

   \frac{\tilde{y}_{j+1} - 2 \tilde{y}_{j} + \tilde{y}_{j-1}}{\Delta x^{2}} = - \frac{\sigma^{2}}{c^{2}} \tilde{y}_{j}
   \qquad (2 \leq j \leq N-1).

Together with the two boundary conditions

.. math::

   \tilde{y}_{1} = 0 \qquad
   \tilde{y}_{N} = 0,

we thus have a linear system of :math:`N` algebraic equations and :math:`N` unknowns.
   
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

   \tau \equiv \frac{\Delta x} c

as the sound crossing time of a single cell.

Equation (:eq:`linear-sys`) is a :wiki:`homogeneous linear system
<System_of_linear_equations#Homogeneous_systems>`, meaning that it
has non-trivial solutions :math:`\vu` only when the determinant of
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
\ldots` predicted by the analytic solutions.

Scanning for Eigenfrequencies
-----------------------------

While :numref:`fig-discrim-func` is useful for visualizing
:math:`\Dfunc`, it's not the best way to find
eigenfrequencies. Instead, we can rely on well-established techniques
for isolating and refining roots of monovariate functions.

First, we evaluate a set of :math:`M` values
:math:`\{\Dfunc_{1},\Dfunc_{2},\ldots,\Dfunc_{M}\}`, representing the
discriminant function sampled on the discrete frequency grid
:math:`\{\sigma_{1},\sigma_{2},\ldots,\sigma_{M}\}`. Then, we scan
through these data looking for sign changes between adjacent
discriminant values. If :math:`\Dfunc_{i} \Dfunc_{i+1} < 0`, we know
that a root of the discriminant function must lie in the interval
:math:`(\sigma_{i},\sigma_{i+1})` --- we have *bracketed* a
root. :numref:`fig-discrim-brackets` illustrates the process of
bracket scanning for a frequency grid comprising :math:`M=32` points,
distributed uniformly in :math:`\sigma` across the same range as
plotted in :numref:`fig-discrim-func`. This figure highlights five
brackets containing the five roots identified previously.

.. _fig-discrim-brackets:

.. figure:: fig_discrim_brackets.svg
   :alt: Plot showing the discriminant function versus frequency, with
         root brackets indicated
   :align: center

   Plot of the discriminant values :math:`\{\Dfunc\}` on the discrete
   frequency grid :math:`\{\sigma\}` (distributed uniformly in
   :math:`\sigma`), for the stretched-string BVP with :math:`N=50` and
   :math:`M=32`. The orange-haloed segments highlight adjacent points
   that bracket a root :math:`\Dfunc=0`. (:download:`Source
   <fig_discrim_brackets.py>`)

Once a bracket is established for a given root, it can be narrowed
through a process of iterative refinement until the root is converged
upon. There are a variety of well-known root-finding algorithms that
perform this refinement; the :wiki:`bisection method <Bisection_method>` is conceptually
the simplest, but approaches such as :wiki:`Brent's method <Brent's_method>` can be
much more efficient. For the brackets plotted in
:numref:`fig-discrim-brackets`, :numref:`numerical-eigenfreqs` compares
the eigenfrequencies found using Python's
:external:py:func:`scipy.optimize.brentq` function, against the analytic values
predicted by equation (:eq:`analytic-eigenfreqs`).

.. _numerical-eigenfreqs:

.. csv-table:: Numerical and analytic eigenfrequencies, in units of
   :math:`\pi c/L`, for the stretched-string BVP with
   :math:`N=50`. (:download:`Source <discrim_roots.py>`)
   :widths: 20 40 40
   :align: center
   :file: discrim_roots.csv

Eigenfunction Reconstruction
----------------------------

For each of the eigenfrequencies found, we reconstruct the
corresponding eigenfunction by solving the linear system
(:eq:`linear-sys`). Because :math:`\det(\mS)` is now zero, this system
is guaranteed to have a non-trivial solution. The solution vector
:math:`\vu` resides in the :wiki:`null space <Null_space>` of
:math:`\mS`, and we can use standard numerical techniques (e.g.,
:wiki:`singular value decomposition <Singular_value_decomposition>`)
to evaluate it.  Then, the :math:`j`'th element of :math:`\vu`
corresponds to the eigenfunction sampled at the :math:`j`'th spatial
grid point:

.. math::

   (\vu)_{j} = \tilde{y}_{j} \equiv \tilde{y}(x_{j})

.. _fig-eigenfuncs:

.. figure:: fig_eigenfuncs.svg
   :alt: Plot showing eigenfunctions for the first three modes
   :align: center

   Plot of the eigenfunctions :math:`\tilde{y}` as a function of
   spatial coordinate :math:`x`, for the first three modes of the
   stretched-string BVP with :math:`N=50`. The discrete points show
   the numerical functions, and the solid lines the corresponding
   analytic functions. In all cases, the eigenfunctions have been
   normalized to have a maximum :math:`|\tilde{y}|` of
   unity. (:download:`Source <fig_eigenfuncs.py>`)

:numref:`fig-eigenfuncs` plots the eigenfunctions found in this way
for the first three modes (:math:`n=1,2,3`) given in
:numref:`numerical-eigenfreqs`. Also shown are the corresponding
analytic solutions given by equation (:eq:`analytic-eigenfuncs`). The
agreement between the two is good.
