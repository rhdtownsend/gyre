.. _osc-dimless-form:

Dimensionless Formulation
=========================

To improve numerical stability, :program:`gyre` solves the
:ref:`separated equations <osc-sep-eqns>` by recasting them into a
dimensionless form that traces its roots back to
:ads_citet:`dziembowski:1971`.

.. _osc-dimless-vars:

Variables
---------

The independent variable is the fractional radius :math:`x \equiv r/R`
and the dependent variables :math:`\{y_{1},y_{2},\ldots,y_{6}\}` are

.. math::
   :label: e:dimless

   \begin{align}
   y_{1} &= x^{2 - \ell}\, \frac{\txir}{r}, \\
   y_{2} &= x^{2-\ell}\, \frac{\tP'}{\rho g r}, \\
   y_{3} &= x^{2-\ell}\, \frac{\tPhi'}{gr}, \\
   y_{4} &= x^{2-\ell}\, \frac{1}{g} \deriv{\tPhi'}{r}, \\
   y_{5} &= x^{2-\ell}\, \frac{\delta \tS}{c_{p}}, \\
   y_{6} &= x^{-1-\ell}\, \frac{\delta \tLrad}{L}.
   \end{align}

.. _osc-dimless-eqns:

Oscillation Equations
---------------------

The dimensionless oscillation equations are

.. math::

   \begin{align}
   x \deriv{y_{1}}{x} &=
   \left( \frac{V}{\Gammi} - 1 - \ell \right) y_{1} +
   \left( \frac{\ell(\ell+1)}{c_{1} \omegac^{2}} - \alphagam \frac{V}{\Gammi} \right) y_{2} +
   \alphagrv \frac{\ell(\ell+1)}{c_{1} \omegac^{2}} y_{3} +
   \delta \, y_{5}, \\
   %
   x \deriv{y_{2}}{x} &=
   \left( c_{1} \omegac^{2} - \fpigam \As \right) y_{1} +
   \left( 3 - U + \As - \ell \right) y_{2} -
   \alphagrv y_{4} +
   \delta \, y_{5}, \\
   %
   x \deriv{y_{3}}{x} &= 
   \alphagrv \left( 3 - U - \ell \right) y_{3} +
   \alphagrv y_{4} \\
   %
   x \deriv{y_{4}}{x} &=
   \alphagrv \As U y_{1} +
   \alphagrv \frac{V}{\Gammi} U y_{2} +
   \alphagrv \ell(\ell+1) y_{3} -
   \alphagrv (U + \ell - 2) y_{4}
   - \alphagrv \delta \, U y_{5}, \\
   %
   x \deriv{y_{5}}{x} &= 
   \frac{V}{\frht} \left[ \nabad (U - c_{1}\omegac^{2}) - 4 (\nabad - \nabla) + \ckapad V \nabla + \cdif \right] y_{1} + \mbox{} \\
   &
   \frac{V}{\frht} \left[ \frac{\ell(\ell+1)}{c_{1} \omegac^{2}} (\nabad - \nabla) - \ckapad V \nabla - \cdif \right] y_{2} + \mbox{} \\
   &
   \alphagrv \frac{V}{\frht} \left[ \frac{\ell(\ell+1)}{c_{1} \omegac^{2}} (\nabad - \nabla) \right] y_{3} +
   \alphagrv \frac{V \nabad}{\frht} y_{4} + \mbox{} \\
   &
   \left[ \frac{V \nabla}{\frht} (4 \frht - \ckapS) + \dfrht + 2 - \ell \right] y_{5} -
   \frac{V \nabla}{\frht \crad} y_{6} \\
   %
   x \deriv{y_{6}}{x} &=
   \left[ \alphahfl \ell(\ell+1) \left( \frac{\nabad}{\nabla} - 1 \right) \crad - V \cepsad - \alphaegv c_{egv} \nabad V \right] y_{1} + \mbox{} \\
   &
   \left[ V \cepsad - \ell(\ell+1) \crad \left( \alphahfl \frac{\nabad}{\nabla} - \frac{3 + \dcrad}{c_{1}\omegac^{2}} \right) + \alphaegv c_{egv} \nabad V \right] y_{2} + \mbox{} \\
   &
   \alphagrv \left[ \ell(\ell+1) \crad \frac{3 + \dcrad}{c_{1}\omegac^{2}} \right] y_{3} + \mbox{} \\
   &
   \left[ \cepsS - \alphahfl \frac{\ell(\ell+1)\crad}{\nabla V} + \ii \alphathm \omegac \cthk + \alphaegv c_{egv} \right] y_{5} -
   \left[ 1 + \ell \right] y_{6}.
   \end{align}

These equations are derived from the separated equations, but with the
insertion of 'switch' terms (denoted :math:`\alpha`) that allow
certain pieces of physics to be altered. See the
:ref:`osc-physics-switches` section for more details

For non-radial adiabatic calculations, the last two equations above
are set aside and the :math:`y_{5}` terms dropped from the first four
equations. For radial adiabatic calculations with
:nml_n:`reduce_order`\ =\ :nml_v:`.TRUE.` (see the :ref:`osc-params`
section), the last four equations are set aside and the first two
replaced by

.. math::

   \begin{align}
   x \deriv{y_{1}}{x} &=
   \left( \frac{V}{\Gammi} - 1 \right) y_{1} - \frac{V}{\Gamma_{1}} y_{2}, \\
   %
   x \deriv{y_{2}}{x} &= 
   \left( c_{1} \omega^{2} + U - \As \right) y_{1} + \left( 3 - U + \As \right) y_{2}.
   \end{align}

.. _osc-dimless-bc:

Boundary Conditions
-------------------

Inner Boundary
^^^^^^^^^^^^^^

When :nml_n:`inner_bound`\ =\ :nml_v:`'REGULAR'`, GYRE applies
regularity-enforcing conditions at the inner boundary:

.. math::

   \begin{align}
   c_{1} \omega^{2} y_{1} - \ell y_{2} - \alphagrv \ell y_{3} &= 0, \\
   \alphagrv \ell y_{3} - (2\alphagrv - 1) y_{4} &= 0, \\
   y_{5} &= 0.
   \end{align}

When :nml_n:`inner_bound`\ =\ :nml_v:`'ZERO_R'`, the first and second
conditions are replaced with zero radial displacement conditions,

.. math::
   
   \begin{align}
   y_{1} &= 0, \\
   y_{4} &= 0.
   \end{align}

Likewise, when :nml_n:`inner_bound`\ =\ :nml_v:`'ZERO_H'`, the first and
second conditions are replaced with zero horizontal displacement
conditions,

.. math::

   \begin{align}
   y_{2} - y_{3} &= 0, \\
   y_{4} &= 0.
   \end{align}

Outer Boundary
^^^^^^^^^^^^^^

When :nml_n:`outer_bound`\ =\ :nml_v:`'VACUUM'`, GYRE applies vacuum surface
pressure conditions at the outer boundary:

.. math::
   :label: e:outer-bc

   \begin{align}
   y_{1} - y_{2} &= 0 \\
   \alphagrv U y_{1} + (\alphagrv \ell + 1) y_{3} + \alphagrv y_{4} &= 0 \\
   (2 - 4\nabad V) y_{1} + 4 \nabad V y_{2} + 4 \frht y_{5} - y_{6} &= 0
   \end{align}
   
When :nml_n:`outer_bound`\ =\ :nml_v:`'DZIEM'`, the first condition is
replaced by the :ads_citet:`dziembowski:1971` outer mechanical
boundary condition,

.. math::

   \left\{ 1 + V^{-1} \left[ \frac{\ell(\ell+1)}{c_{1} \omega^{2}} - 4 - c_{1} \omega^{2} \right] \right\} y_{1} -
   y_{2} = 0.
   
When :nml_n:`outer_bound`\ =\ :nml_v:`'UNNO'`\ \|\ :nml_v:`'JCD'`, the
first condition is replaced by the (possibly-leaky) outer mechanical
boundary conditions described by :ads_citet:`unno:1989` and
:ads_citet:`christensen-dalsgaard:2008`, respectively. When
:nml_n:`outer_bound`\ =\ :nml_v:`'ISOTHERMAL'`, the first condition is
replaced by a (possibly-leaky) outer mechanical boundary condition
derived from a local dispersion analysis of an isothermal atmosphere.

Finally, when :nml_n:`outer_bound`\ =\ :nml_v:`'GAMMA'`, the first
condition is replaced by the outer mechanical boundary condition
described by :ads_citet:`ong:2020`.

.. _osc-dimless-jump:

Jump Conditions
---------------
  
Across density discontinuities, GYRE enforces conservation of mass,
momentum and energy by applying the jump conditions

.. math::
   
   \begin{align}
   U^{+} y_{2}^{+} - U^{-} y_{2}^{-} &= y_{1} (U^{+} - U^{-}) \\
   y_{4}^{+} - y_{4}^{-} &= -y_{1} (U^{+} - U^{-}) \\
   y_{5}^{+} - y_{5}^{-} &= - V^{+} \nabad^{+} (y_{2}^{+} - y_{1}) +
   V^{-} \nabad^{-} (y_{2}^{-} - y_{1})
   \end{align}

Here, + (-) superscripts indicate quantities evaluated on the inner
(outer) side of the discontinuity. :math:`y_{1}`, :math:`y_{3}` and
:math:`y_{6}` remain continuous across discontinuities, and therefore
don't need these superscripts.

.. _osc-struct-coeffs:

Structure Coefficients
----------------------

The various stellar structure coefficients appearing in the
dimensionless oscillation equations are defined as follows:

.. math::

   \begin{gather}
   V = -\deriv{\ln P}{\ln r} \qquad
   V_{2} = x^{-2} V \qquad
   \As = \frac{1}{\Gamma_{1}} \deriv{\ln P}{\ln r} - \deriv{\ln \rho}{\ln r} \qquad
   U = \deriv{\ln M_{r}}{\ln r} \\
   %
   c_1 = \frac{r^{3}}{R^{3}} \frac{M}{M_{r}} \qquad
   \fpigam =
   \begin{cases}
   \alphapi & \As > 0, x < x_{\rm atm} \\
   \alphagam & \As > 0, x > x_{\rm atm} \\
   1 & \text{otherwise}
   \end{cases}\\
   %
   \nabla = \deriv{\ln T}{\ln P} \qquad
   \clum = x^{-3} \frac{\Lrad+\Lcon}{L} \qquad
   \crad = x^{-3} \frac{\Lrad}{L} \qquad
   \dcrad = \deriv{\ln \crad}{\ln r} \\
   %
   \frht = 1 - \alpharht \frac{\ii \omega \cthn}{4} \qquad
   \dfrht = - \alpharht \frac{\ii \omega \cthn \dcthn}{4 \frht} \\
   %
   \ckapad = \frac{\alphakar \kaprho}{\Gamma_{1}} + \nabad \alphakat \kapT \qquad
   \ckapS = - \upsT \alphakar \kaprho + \alphakat \kapT \\
   %
   \ceps = x^{-3} \frac{4\pi r^{3} \rho \epsnuc}{L} \qquad
   \cepsad = \ceps \epsnucad \qquad
   \cepsS = \ceps \epsnucS \\
   %
   \cdif = - 4 \nabad V \nabla + \nabad \left(V + \deriv{\ln \nabad}{\ln x} \right) \\
   %
   \cthn = \frac{\cP}{a c \kappa T^{3}} \sqrt{\frac{GM}{R^{3}}} \qquad
   \dcthn = \deriv{\ln \cthn}{\ln r} \\
   %
   \cthk = x^{-3} \frac{4\pi r^{3} \cP T \rho}{L} \sqrt{\frac{GM}{R^{3}}}
   \cegv = x^{-3} \frac{4\pi r^{3} \rho \epsgrav}{L}
   \end{gather}

.. _osc-physics-switches:

Physics Switches
----------------

GYRE offers the capability to adjust the oscillation equations through
a number of physics switches, controlled by parameters in the
:nml_g:`osc` namelist group. The table below summarizes the mapping
between the switches appearing in the expressions above, and the
corresponding namelist parameters.

.. list-table::
   :widths: 20 20 60
   :header-rows: 1

   * - Symbol
     - Parameter
     - Description
   * - :math:`\alphagrv`
     - :nml_n:`alpha_grv`
     - Scaling factor for gravitational potential perturbations. Set to 1
       for normal behavior, and to 0 for the :ads_citet:`cowling:1941`
       approximation
   * - :math:`\alphathm`
     - :nml_n:`alpha_thm`
     - Scaling factor for local thermal timescale. Set to 1 for normal
       behavior, to 0 for the non-adiabatic reversible (NAR) approximation
       (see :ads_citealp:`glatzel:1990`), and to a large value to approach
       the adiabatic limit
   * - :math:`\alphahfl`
     - :nml_n:`alpha_hfl`
     - Scaling factor for horizontal flux perturbations. Set to 1 for
       normal behavior, and to 0 for the non-adiabatic radial flux (NARF)
       approximation (see :ads_citealp:`townsend:2003b`)
   * - :math:`\alphagam`
     - :nml_n:`alpha_gam`
     - Scaling factor for g-mode isolation. Set to 1 for normal behavior,
       and to 0 to isolate g modes as described by :ads_citet:`ong:2020`
   * - :math:`\alphapi`
     - :nml_n:`alpha_pi`
     - Scaling factor for p-mode isolation. Set to 1 for normal behavior,
       and to 0 to isolate p modes as described by :ads_citet:`ong:2020`
   * - :math:`\alphakar`
     - :nml_n:`alpha_kar`
     - Scaling factor for opacity density partial derivative. Set to 1 for normal
       behavior, and to 0 to suppress the density part of the :math:`\kappa` mechanism
   * - :math:`\alphakat`
     - :nml_n:`alpha_kat`
     - Scaling factor for opacity temperature partial derivative. Set to 1 for normal
       behavior, and to 0 to suppress the temperature part of the :math:`\kappa` mechanism
   * - :math:`\alpharht`
     - :nml_n:`alpha_rht`
     - Scaling factor for time-dependent term in the radiative heat
       equation (see :ads_citealp:`unno:1966`). Set to 1 to include this
       term (Unno calls this the Eddington approximation), and to 0 to
       ignore the term
   * - :math:`\alphaegv`
     - :nml_n:`alpha_egv`
     - Scaling factor for the gravitational heating rate terms.
       Set to 1 to include these terms (default), and to 0 to suppress it.
