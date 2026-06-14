.. _osc-tidal:

Tidal Effects
=============

To simulate the effects of tidal forcing by a companion, the
:program:`gyre_tides` frontend solves a modified form of the
linearized momentum equation (:eq:`e:osc-lin-mom`), namely

.. math::

   \rho \pderiv{\vv'}{t} = - \nabla P' - \rho' \nabla P - \rho \nabla \Phi' - \rho \nabla \PhiT.

The final term on the right-hand side represents the external force
arising from the tidal gravitational potential :math:`\PhiT`.

Tidal Potential
---------------

The tidal potential can be expressed via the superposition

.. math::
   :label: e:tidal-pot

   \PhiT = \sum_{\ell=2}^{\infty} \sum_{m=-\ell}^{\ell} \sum_{k=-\infty}^{\infty} \PhiTlmk.

of partial tidal potentials defined by

.. math::
   :label: e:partial-tidal-pot

   \PhiTlmk \equiv
   - \epsT \,
   \frac{G\Mstar}{\Rstar} \,
   \cbar_{\ell,m,k}
   \left( \frac{r}{\Rstar} \right)^{\ell} Y^{m}_{\ell}(\theta, \phi) \,
   \exp(- \ii k \Oorb t).

(the summation over :math:`\ell` and :math:`m` comes from a multipolar
space expansion of the potential, and the summation over :math:`k`
from a Fourier time expansion). Here,

.. math::

   \epsT = \left( \frac{\Rstar}{a} \right)^{3} q_{2} = \frac{\Oorb^{2} \Rstar^{3}}{G\Mstar} \frac{q_{2}}{1+q_{2}}

quantifies the overall strength of the tidal forcing, in terms of the
companion's mass :math:`q_{2} \Mstar`, semi-major axis :math:`a` and orbital
angular frequency :math:`\Oorb`. These expressions, and the definition
of the tidal expansion coefficients :math:`\cbar_{\ell,m,k}`, are presented in
greater detail in :ads_citet:`sun:2023`.

Tidal Equations
---------------

Because the tidal potential (:eq:`e:tidal-pot`) superposes many
different spherical harmonics, the solution forms given in equations
(:eq:`e:osc-sol-form`) must be replaced by

.. math::
   :label: e:osc-sol-form-tide

   \begin{aligned}
   f'(r,\theta,\phi;t) &=
   \sum_{\ell,m,k} \tflmk'(r) \, Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii k \Oorb t) \\
   \va'(r,\theta,\phi;t) &=
   \sum_{\ell,m,k} \left[ \tarlmk'(r) \, \ver + \tahlmk'(r) \, r \, \nablah \right] Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t),
   \end{aligned}

and those in equation (:eq:`e:osc-sol-form-vxi`) by

.. math::
   :label: e:osc-sol-form-vxi-tide

   \vxi(r,\theta,\phi;t) =
   \sum_{\ell,m,k} \left[ \txirlmk(r) \, \ver + \txihlmk(r) \, r \, \nablah \right] Y^{m}_{\ell}(\theta,\phi) \, \exp(-\ii \sigma t)
	   
(the notation for the sums has been abbreviated). Substituting these
solution forms into the :ref:`linearized equations <osc-linear-eqns>`
(with the additional :math:`\rho \nabla \PhiT` term, as above), and
taking advantage of the orthonormality of the spherical harmonics,
leads to a fully separated set of differential equations for each
combination of :math:`\ell`, :math:`m` and :math:`k`. A given set
resembles the regular :ref:`separated equations <osc-sep-eqns>`, but
with the radial momentum equation (:eq:`e:osc-sep-mom-r`) replaced by

.. math::

   -\sigmac^{2} \rho \txir = - \deriv{\tP'}{r} - \trho' \deriv{\Phi}{r} - \rho \deriv{\tPhi'}{r} - \rho \deriv{\tPhiT}{r}

and the horizontal momentum equation (:eq:`e:osc-sep-mom-h`) by 

.. math::

   -\sigmac^{2} \rho \txih = - \frac{1}{r} \left( - \tP' - \rho \tPhi' - \rho \tPhiT \right).

(in the interests of brevity, the :math:`\ell,m,k` are now
dropped). Here,

.. math::

   \tPhiT =  - \epsT \, \frac{G\Mstar}{\Rstar} \, \cbar_{\ell,m,k} \left( \frac{r}{\Rstar} \right)^{\ell}

describes the radial dependence of the partial tidal potential defined
by eqn. (:eq:`e:partial-tidal-pot`), while

.. math::

   \sigmac = k \Oorb - m \Orot

is the forcing frequency associated with the potential, corrected for
the :ref:`Doppler shift <osc-rot-doppler>` due to rotation.

In the :ref:`dimensionless oscillation equations <osc-dimless-eqns>`,
the tidal modifications consist of additional inhomogeneous terms on
the right-hand sides of the differential equations for :math:`y_1`,
:math:`y_2`, :math:`y_5` and :math:`y_6`, as follows:

.. math::

   \begin{aligned}
   x \deriv{y_{1}}{x} &= [ \ldots ] + \frac{\lambda}{c_{1} \omegac^{2}} \yT \\
   x \deriv{y_{2}}{x} &= [ \ldots ] - \ell \yT \\
   x \deriv{y_{5}}{x} &= [ \ldots ] + \frac{V}{\frht} \left[ \frac{\lambda}{c_{1} \omegac^{2}} (\nabad - \nabla) + \ell \nabad \right] \yT \\
   x \deriv{y_{6}}{x} &= [ \ldots ] + \left[ \lambda \crad \frac{3 + \dcrad}{c_{1}\omegac^{2}} \right] \yT
   \end{aligned}

where :math:`[\ldots]` represents the right-hand side in the absence
of tidal forcing, and

.. math::

   \yT \equiv x^{2-\ell}\, \frac{\tPhiT}{gr} = - \epsT \, c_{1} \cbar_{\ell,m,k}.

(The differential equations for :math:`y_{3}` and :math:`y_{4}` remain
unchanged). Moreover, the first of the regularity-enforcing equations
(:eq:`e:inner-bound`) becomes

 .. math::

   c_{1} \omegac^{2} y_{1} - \ell y_{2} - \alphagrv \ell y_{3} = \ell \yT.

In these expressions,

.. math::

   \omegac = \alphafrq \sigmac \sqrt{\frac{\Rstar^{3}}{G\Mstar}}

is the dimensionless tidal forcing frequency in the co-rotating
reference frame. The :math:`\alphafrq` term is included here to allow
the forcing frequency to be adjusted independently of :math:`m` and
:math:`k` (see the :nml:option:`alpha_frq <tide.alpha_frq>`
option). All occurrences of :math:`\omega` in the other dimensionless
equations are likewise replaced by :math:`\omegac`.
