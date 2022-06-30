.. _tide-params:

Tide Parameters
===============

The :nml_g:`tide` namelist group defines tidal parameters, as follows:

:nml_n:`Phi_T_thresh` (default :nml_v:`0.`)
  Threshold on tidal potential :math:`Phi_{\rm T}` for a component tide to contribute to tidal forcing

:nml_n:`omega_c_thresh` (default :nml_v:`0.`) Threshold on
  dimensionless co-rotatign frequency :math:`\omega_{\rm c}` for a
  component tide to be treated as dynamic (rather than static)

:nml_n:`alpha_frq` (default :nml_v:`1.`)
  Scaling parameter :math:`\alphafrq` for tidal forcing frequency

:nml_n:`l_min` (default :nml_v:`2`)
  Minimum harmonic degree :math:`\ell` in spatial expansion of tidal potential

:nml_n:`l_max` (default :nml_v:`2`)
  Maximum harmonic degree :math:`\ell` in spatial expansion of tidal potential

:nml_n:`m_min` (default :nml_v:`-HUGE`)
  Minimum azimuthal order :math:`m` in spatial expansion of tidal potential

:nml_n:`m_max` (default :nml_v:`HUGE`)
  Maximum azimuthal order :math:`m` in spatial expansion of tidal potential

:nml_n:`k_min` (default :nml_v:`-10`)
  Minimum orbital harmonic :math:`k` in temporal expansion of tidal potential

:nml_n:`k_max` (default :nml_v:`10`)
  Maximum orbital harmonic :math:`k` in temporal expansion of tidal potential

:nml_n:`tag`
  Tag for controlling selection of other parameters

.. rubric:: Footnotes

.. [#only-D] This option is available only for stellar models with :ref:`D capability <model-caps>`
