.. _orbit-group:

.. nml:group:: orbit

Orbit Namelist Group
====================

The :nml:group:`orbit` namelist group defines orbital parameters. The
input file can contain one or more, but only the last (tag-matching)
one is used. The following options are available:

.. nml:option:: Omega_orb
   :type: real
   :default: 1

   Orbital angular frequency

.. nml:option:: Omega_orb_units
   :type: string
   :default: 'NONE'

   Units of :nml:option:`Omega_orb`; one of

   - :nml:value:`'NONE'` : Dimensionless angular frequency
   - :nml:value:`'HZ'` : Linear frequency in Hz\ [#only-D]_
   - :nml:value:`'UHZ'` : Linear frequency in :math:`\mu`\ Hz\ [#only-D]_
   - :nml:value:`'RAD_PER_SEC'` : Angular frequency in radians per second\ [#only-D]_
   - :nml:value:`'CYC_PER_DAY'` : Linear frequency in cycles per day\ [#only-D]_
   - :nml:value:`'CRITICAL'` : Fraction of the Roche critical rate\ [#only-D]_

.. nml:option:: q
   :type: real
   :default: 1

   Ratio of secondary mass to primary mass

.. nml:option:: e
   :type: real
   :default: 0

   Orbital eccentricity

.. nml:option:: tag_list
   :type: string
   :default: ''

   Comma-separated list of :nml:option:`tag <mode.tag>` values to match;
   matches all if left blank

.. rubric:: Footnotes

.. [#only-D] This option is available only for stellar models with :ref:`D capability <model-caps>`
