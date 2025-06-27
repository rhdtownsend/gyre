.. _osc-rot:

Rotation Effects
================

The oscillation equations presented in the preceding sections are
formulated for a non-rotating star. The corresponding equations for a
rotating star are significantly more complicated, and a complete
treatment of rotation lies beyond the scope of GYRE. However, GYRE can
include two important effects arising from rotation.

.. _osc-rot-doppler:

Doppler Shift
-------------

A lowest-order effect of rotation arises in the Doppler shift from
transforming between the inertial reference frame and the local
co-rotating reference frame. To incorporate this effect in the
:ref:`separated equations <osc-sep-eqns>`, all instances of the
inertial-frame frequency :math:`\sigma` are replaced by the
co-rotating frequency

.. math::
   :label: e:sigmac

   \sigmac \equiv \sigma - m \Orot,

where :math:`m` is the azimuthal order of the mode and :math:`\Orot`
is the rotation angular frequency. GYRE assumes shellular rotation
(see, e.g., :ads_citealp:`meynet:1997`), and so the latter can in
principle be a function of radial coordinate :math:`r`. The
corresponding modifications to the :ref:`dimensionless formulation
<osc-dimless-form>` involve replacing the dimensionless inertial-frame
frequency :math:`\omega` with the dimensionless co-rotating frequency

.. math::

   \omegac \equiv \omega - m \Orot \sqrt{\frac{\Rstar^{3}}{G\Mstar}}.

.. _osc-rot-coriolis-p:

Perturbative Coriolis Force Treatment
-------------------------------------

Another lowest-order effect of rotation arises from the Coriolis
force. For slow rotation, this effect can be determined through a
perturbation expansion technique (see, e.g., section 19.2 of
:ads_citealp:`unno:1989`). To first order in :math:`\Orot`, the
frequency of a mode is shifted by the amount

.. math::

   \Delta \sigma = m \int_{0}^{\Rstar} \Orot \, \deriv{\beta}{r} \diff{r},

where the rotation splitting kernel is

.. math::

   \deriv{\beta}{r} =
   \frac{\left\{ \txir^{2} + [\ell(\ell+1) - 1] \txih^{2} - 2 \txir \txih \right\} \rho r^{2}}
   {\int_{0}^{\Rstar} \left\{ \txir^{2} + \ell(\ell+1) \txih^{2} \right\} \rho r^{2} \diff{r}}

In this latter expression, the eigenfunctions :math:`\txir` and
:math:`\txih` are evaluated from solutions to the oscillation
equations without rotation. Therefore, the expression above for
:math:`\Delta \sigma` can be applied as a post-calculation correction
to non-rotating eigenfrequencies.

.. _osc-rot-coriolis-np:

Non-Perturbative Coriolis Force Treatment
-----------------------------------------

The perturbation expansion technique above breaks down when
:math:`\Orot/\sigmac \gtrsim 1`. To deal with such cases, the
:program:`gyre` frontend [#gyre-tides]_ can incorporate a
non-perturbative treatment of the Coriolis force based on the
'traditional approximation of rotation' (TAR). The TAR was first
introduced by Eckart (1960; `Hydrodynamics of Oceans and Atmospheres`)
and has since been used extensively within the pulsation community
(see, e.g., :ads_citealp:`bildsten:1996`; :ads_citealp:`lee:1997`;
:ads_citealp:`townsend:2003a`; :ads_citealp:`bouabid:2013`;
:ads_citealp:`townsend:2020`).

Within the TAR, the solution forms given in
equation (:eq:`e:osc-sol-forms`) are replaced by

.. math::
   :label: e:osc-sol-forms-hough

   \begin{aligned}
   \xir(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \txir(r) \, \houghr(\theta) \, \exp(\ii m \phi -\ii \sigma t) \right], \\
   \xit(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \txih(r) \, \frac{\hought(\theta)}{\sin\theta} \, \exp(\ii m \phi -\ii \sigma t) \right], \\
   \xip(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \txih(r) \, \frac{\houghp(\theta)}{\ii \sin\theta} \, \exp(\ii m \phi -\ii \sigma t) \right], \\
   f'(r,\theta,\phi;t) &= \operatorname{Re} \left[ \sqrt{4\pi} \, \tf'(r) \, \houghr(\theta) \, \exp(\ii m \phi -\ii \sigma t) \right]
   \end{aligned}

Here, the Hough functions :math:`\houghr`, :math:`\hought` and
:math:`\houghp` are the eigenfunctions obtained by solving Laplace's
tidal equations (TEs), a second-order system of differential equations
and boundary conditions in the polar (:math:`\theta`) coordinate (see
:ads_citealt:`townsend:2020`). Together with the associated eigenvalue
:math:`\lambda`, they depend on the harmonic degree :math:`\ell`\
[#harmonic-deg]_ and azimuthal order :math:`m`, and the spin parameter

.. math::

   q \equiv \frac{2 \Orot}{\sigmac}.

.. _osc-rot-solfam:

Solution Families
^^^^^^^^^^^^^^^^^

Solutions to the TEs can be grouped into two families based on the
behavior of the eigenfunctions and eigenvalue in the limit :math:`\Orot
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
   \text{as } \Orot \rightarrow 0.

and :math:`\lambda \rightarrow \ell(\ell+1)`. With these expressions,
the solution forms (:eq:`e:osc-sol-forms-hough`) reduce to those given
in equation (:eq:`e:osc-sol-forms`).

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
   \text{as } \Orot \rightarrow 0.

and :math:`\lambda \rightarrow 0`. Moreover, Rossby-mode
eigenfrequencies also show the limiting behavior

.. math::
   :label: e:ross-freq

   \sigmac = \frac{2 m \Orot}{\ell(\ell+1)}
   \quad
   \text{as } \Orot \rightarrow 0,

which is independent of the stellar structure.

Implementing the TAR
^^^^^^^^^^^^^^^^^^^^

To implement the TAR in the :ref:`separated equations
<osc-sep-eqns>` and :ref:`boundary conditions <osc-bound-conds>`,
all instances of the term :math:`\ell(\ell+1)` are replaced by the TE
eigenvalue :math:`\lambda`. Then, all instances of the harmonic degree
:math:`\ell` are replaced by :math:`\elle`, an effective harmonic
degree found by solving

.. math::

   \elle(\elle+1) = \lambda.

Similar steps are taken in the :ref:`dimensionless formulation
<osc-dimless-form>`, but in the definitions of the dependent variables
:math:`\{y_{1},y_{2},\ldots,y_{6}\}`, :math:`\ell` is replaced by
:math:`\elli`, the effective harmonic degree evaluated at the inner
boundary.

.. rubric:: Footnotes

.. [#gyre-tides] Currently the TAR cannot be used with the
		 :program:`gyre_tides` frontend, because it doesn't play well with
		 forcing by the tidal potential :math:`\PhiT`.

.. [#harmonic-deg] The harmonic degree isn't formally a 'good' quantum
                   number in the TAR; however, it can still be used to
                   identify Hough functions by considering their
                   behavior in the limit :math:`\Orot \rightarrow 0`,
                   as given in eqns. (:eq:`e:hough-lim-ga`) and
                   (:eq:`e:hough-lim-ross`).
