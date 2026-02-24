.. _tide-group:

.. nml:group:: tide

Tide Namelist Group
===================

The :nml:group:`tide` namelist group defines tidal indices and related
quantities. The following options are available:

.. nml:option:: y_T_thresh_abs
   :type: real
   :default: 0

   Absolute threshold on dimensionless tidal potential :math:`y_{\rm
   T}` in order for a partial tide to contribute

.. nml:option:: y_T_thresh_rel
   :type: real
   :default: 0

   Relative threshold on dimensionless tidal potential :math:`y_{\rm
   T}` in order for a partial tide to contribute

.. nml:option:: omega_c_thresh
   :type: real
   :default: 0

   Threshold on dimensionless co-rotating frequency :math:`\omegac` in
   order for a partial tide to be treated as dynamic rather than
   static

.. nml:option:: alpha_frq
   :type: real
   :default: 1

   Scaling factor for the dimensionless co-rotating frequerncy
   :math:`\omegac` (see eqn. :eq:`e:omegac-force`)

.. nml:option:: l_min
   :type: integer
   :default: 2

   Minimum harmonic degree :math:`\ell` in spatial expansion of tidal
   potential

.. nml:option:: l_max
   :type: integer
   :default: 2

   Maximum harmonic degree :math:`\ell` in spatial expansion of tidal
   potential

.. nml:option:: m_min
   :type: integer
   :default: -HUGE

   Minimum azimuthal order :math:`m` in spatial expansion of tidal
   potential

.. nml:option:: m_max
   :type: integer
   :default: HUGE

   Maximum azimuthal order :math:`m` in spatial expansion of tidal
   potential

.. nml:option:: k_min
   :type: integer
   :default: -10

   Minimum orbital harmonic :math:`k` in temporal expansion of tidal
   potential

.. nml:option:: k_max
   :type: integer
   :default: 10

   Maximum orbital harmonic :math:`k` in temporal expansion of tidal
   potential

.. nml:option:: tag
   :type: string

   Tag for controlling selection of other namelist groups (see
   :ref:`working-with-tags`)
