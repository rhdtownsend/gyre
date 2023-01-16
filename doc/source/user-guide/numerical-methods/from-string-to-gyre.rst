.. _numerical-gyre:

From Stretched String to gyre
=============================

The numerical technique demonstrated in the :ref:`numerical-string`
section provides a powerful analog to how :program:`gyre` solves the
oscillation equations. The full details of :program:`gyre`'s approach
are laid out in :ads_citet:`townsend:2013`; in this section we briefly
summarize it, highlighting similarities and differences with the
stretched-string problem.
	   
Separation
----------

Similar to the stretched-string problem, :program:`gyre`
begins by separating variables in space and time. For the radial
displacement perturbation :math:`\xir`, trial solutions take the
form

.. math::

  \xir(r,\theta,\phi;t) = \operatorname{Re} \left[ \sqrt{4\pi} \, \txir(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t) \right]

(this is taken from the :ref:`osc-sep-eqns` section). In addition to
the same sinusoidal time dependence as in eqn. (:eq:`var-separation`), a
spherical harmonic term :math:`Y^{m}_{\ell}` appears because we are
separating in three (spherical) spatial coordinates rather than one.

Discretization
--------------

As with the stretched-string problem, :program:`gyre` discretizes the
ODE governing :math:`\txir(r)` and related quantities on a spatial
grid :math:`\{x_{1},x_{2},\ldots,x_{N}\}`. However, a couple of
important differences arise at this juncture. First, the oscillation
equations are fourth order (sixth, in the non-adiabatic case). Rather
than employing finite-difference approximations to high-order
differential operators, :program:`gyre` instead decomposes the problem
into a system of coupled first-order equations. This system is written
generically as

.. math::

   x \deriv{\vty}{x} = \mA \, \vty,

where :math:`\vty` is a vector of :math:`\neq` dependent variables, and
:math:`\mA` is a :math:`\neq \times \neq` Jacobian matrix. In the
adiabatic case, :math:`\neq=4`; in the non-adiabatic case,
:math:`\neq=6`.

Second, while the above equation system can be discretized using a
simple finite-difference approximation to the left-hand side,
:program:`gyre` offers more-sophisticated approaches with higher
orders of accuracy. These include the Magnus schemes described in
:ads_citet:`townsend:2013`, and implicit Runge-Kutta schemes mentioned
in :ads_citet:`townsend:2018`. The choice of scheme is set by the
:nml_n:`diff_scheme` parameter of the :nml_g:`num` namelist group. The
discretization leads to difference equations of the form

.. math::

   \vty_{j+1} = \mY_{j+1;j} \, \vty_{j},

relating the dependent variable vector at adjacent grid points. The
:math:`\neq \times \neq` fundamental solution matrix :math:`\mY_{j+1,j}`
is evaluated from the value(s) of :math:`\mA` within the interval
:math:`[x_{j},x_{j+1}]` using the discretization scheme.

There are :math:`N-1` of these sets of difference equations. They are
augmented with the boundary conditions

.. math::

   \subin{\mB} \, \vty_{1} = 0,
   \qquad\qquad
   \subout{\mB} \, \vty_{N} = 0,

where :math:`\subin{\mB}` is a :math:`\nin \times \neq` matrix
representing the :math:`\nin` inner boundary conditions, and
:math:`\subout{\mB}` is a :math:`\nout \times \neq` matrix representing
the outer boundary conditions (note that :math:`\nin + \nout =
\neq`). Together, the difference equations and boundary conditions
comprise a linear system of :math:`\neq\,N` algebraic equations
and :math:`\neq N` unknowns.

Linear System
-------------

The linear system can be written in the same form
(cf. eqn. :eq:`linear-sys`) as with the stretched-string problem
. However, now :math:`\vu` is the vector with components

.. math::

   \vu = 
   \begin{pmatrix}
   \vty_{1} \\
   \vty_{2} \\
   \vdots \\
   \vty_{N-1} \\
   \vty_{N}
  \end{pmatrix}

and the system matrix :math:`\mS` is an :math:`\neq N \times \neq N`
block-staircase matrix with components

.. math::

   \mS = 
   \begin{pmatrix}
   \subin{\mB} & \mz & \cdots & \mz & \mz \\
   -\mY_{2;1} & \mI & \cdots & \mz & \mz \\
   \vdots & \vdots & \ddots & \vdots & \vdots \\
   \mz & \mz & \cdots & -\mY_{N;N-1} & \mI \\
   \mz & \mz & \cdots & \mz & \subout{\mB}
   \end{pmatrix}.

As before, the linear system (:eq:`linear-sys`) has non-trivial
solutions only when the determinant of :math:`\mS` vanishes. Thus,
:program:`gyre` finds eigenvalues of the oscillation equation by solving the
characteristic equation

.. math::

   \Dfunc(\omega) \equiv \det(\mS) = 0,

where the dimensionless frequency

.. math::

   \omega \equiv \sqrt{\frac{R^{3}}{GM}} \, \sigma,

is the product of the star's dynamical timescale and the oscillation
frequency :math:`\sigma`. (Internally, :program:`gyre` works
extensively with such :ref:`dimensionless quantities
<osc-dimless-form>`, as it improves the stability of the numerical
algorithms).

Scanning for Eigenfrequencies
-----------------------------

In the adiabatic case, :program:`gyre` searches for roots of the discriminant
function :math:`\Dfunc` using the same bracketing and refinement
strategies as the stretched-string problem.

In the non-adiabatic case, a complication is that the discriminant
function and the dimensionless frequency are both complex
quantities. Solving the characteristic equation in the complex plane
is computationally challenging because there is no equivalent to
bracketing and refinement. :program:`gyre` implements a couple of different
approaches to the problem, as discussed in the :ref:`non-ad-osc`
section.
