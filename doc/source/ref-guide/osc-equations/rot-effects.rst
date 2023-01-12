.. _osc-rot:

Rotation Effects
================

The oscillation equations and boundary conditions laid out in the
:ref:`osc-dimless-form` section are formulated for a non-rotating
star. Solving the corresponding equations for a rotating star is a
challenging task, and a complete treatment lies beyond the scope of
:program:`gyre`. However, :program:`gyre` does include two important
modifications arising from rotation.

.. _osc-rot-doppler:

Doppler Shift
-------------

The lowest-order effect of rotation appears in the Doppler shift that
arises when transforming between the inertial reference frame and the
local co-rotating reference frame. To incorporate this effect in the
:ref:`oscillation equations <osc-sep-eqns>`, all instances of the
inertial-frame frequency :math:`\sigma` are replaced by the
co-rotating frequency

.. math::

   \sigmac \equiv \sigma - m \Omega,

where :math:`m` is the azimuthal order of the mode and :math:`\Omega`
is the rotation angular frequency. Because GYRE assumes so-called
shellular rotation (see, e.g., :ads_citealp:`meynet:1997`), both
:math:`\Omega` and :math:`\sigmac` are functions of radial coordinate
:math:`r`.

.. _osc-rot-coriolis:

Coriolis Force
--------------

.. _osc-rot-tar:

The Traditional Approximation of Rotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Higher-order effects of rotation arise through the Coriolis force,
which appears in the linearized momentum equation to correct for the
non-inertial nature of the co-rotating reference frame. GYRE
incorporates an approximate treatment of the Coriolis force based on
the `traditional approximation of rotation` (TAR), which was first
introduced by Eckart (1960; `Hydrodynamics of Oceans and Atmospheres`)
and has since then been used extensively within the pulsation
community (see, e.g., :ads_citealp:`bildsten:1996`;
:ads_citealp:`lee:1997`; :ads_citealp:`townsend:2003a`;
:ads_citealp:`bouabid:2013`; :ads_citealp:`townsend:2020`).

Within the TAR, the solution forms given in equation
:eq:`e:osc-sol-forms` are replaced by

.. math::
   :label: e:osc-sol-forms-hough

   \begin{aligned}
   \xir(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \txir(r) \, \houghr(\theta) \, \exp(\ii m \phi -\ii \sigma t) \right], \\
   \xit(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \txih(r) \, \frac{\hought(\theta)}{\sin\theta} \, \exp(\ii m \phi -\ii \sigma t) \right], \\
   \xip(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \txih(r) \, \frac{\houghp(\theta)}{\ii \sin\theta} \, \exp(\ii m \phi -\ii \sigma t) \right], \\
   f'(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \tf'(r) \, \houghr(\theta) \, \exp(\ii m \phi -\ii \sigma t) \right]
   \end{aligned}

(cf. equations 1-3 of :ads_citealp:`townsend:2020`). Here, the Hough
functions :math:`\houghr`, :math:`\hought` and :math:`\houghp` are the
eigenfunctions obtained by solving Laplace's tidal equations (TEs), a
second-order system of differential equations and boundary conditions
in the polar (:math:`\theta`) coordinate. Together with the associated
eigenvalue :math:`\lambda`, depend on the harmonic degree
:math:`\ell`\ [#harmonic-deg]_ and azimuthal order :math:`m`, and the
spin parameter

.. math::

   q \equiv \frac{2 \Omega}{\sigmac}.

Solution Families
~~~~~~~~~~~~~~~~~

Solutions to the TEs can be grouped into two families based on the
behavior of the eigenfunctions and eigenvalue in the limit :math:`\Omega
\rightarrow 0`. For the gravito-acoustic family,

.. math::
   :label: e:hough-lim-ga

   \left.
   \begin{aligned}
   \houghr(\theta) \ \rightarrow & \ Y^{m}_{\ell}(\theta,0) \\
   \hought(\theta) \ \rightarrow & \ \sin\theta \pderiv{}{\theta} Y^{m}_{\ell}(\theta,0) \\
   \houghp(\theta) \ \rightarrow & \ - m Y^{m}_{\ell}(\theta,0)
   \end{aligned}
   \right\}
   \quad
   \text{as } \Omega \rightarrow 0.

and :math:`\lambda \rightarrow \ell(\ell+1)`. With these expressions,
the solution forms in equation :eq:`e:osc-sol-forms-hough` reduce to those
given in equation :eq:`e:osc-sol-forms`.

Conversely, for the Rossby family

.. math::
   :label: e:hough-lim-ross

   \left.
   \begin{aligned}
   \houghr(\theta) \ \rightarrow & \ 0 \\
   \hought(\theta) \ \rightarrow & \ m Y^{m}_{\ell}(\theta,0) \\
   \houghp(\theta) \ \rightarrow & \ - \sin\theta \pderiv{}{\theta} Y^{m}_{\ell}(\theta,0)
   \end{aligned}
   \right\}
   \quad
   \text{as } \Omega \rightarrow 0.

and :math:`\lambda \rightarrow 0`. Moreover, Rossby-mode
eigenfrequencies also show the limiting behavior

.. math::
   :label: e:ross-freq

   \sigmac = \frac{2 m \Omega}{\ell(\ell+1)}
   \quad
   \text{as } \Omega \rightarrow 0,

which is independent of the stellar structure.

Incorporating the TAR
~~~~~~~~~~~~~~~~~~~~~

To incorporate the TAR in the oscillation equations, all instances of
the term :math:`\ell(\ell+1)` are replaced by the TE eigenvalue
:math:`\lambda`. Then, all instances of the harmonic degree
:math:`\ell` are replaced by :math:`\elli`, an effective harmonic
degree found by solving

.. math::

   \elli(\elli+1) = \lambda

`at the inner boundary` (remember, because :math:`\sigmac` is a
function of radial coordinate, so too are :math:`q` and
:math:`\lambda`).

.. rubric:: Footnotes

.. [#harmonic-deg] The harmonic degree isn't formally a 'good' quantum
                   number in the TAR; however, it can still be used to
                   identify Hough functions by considering their
                   behavior in the limit :math:`\Omega \rightarrow 0`,
                   as given in eqns. :eq:`e:hough-lim-ga` and
                   :eq:`e:hough-lim-ross`.
