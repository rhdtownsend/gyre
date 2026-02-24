.. _rot-group:

.. nml:group:: rot

Rot Namelist Group
==================

The :nml:group:`rot` namelist group governs the treatment of
rotation. The input file can contain one or more, but only the last
(tag-matching) one is used. The following options are available:

.. nml:option:: Omega_rot_source
   :type: string
   :default: 'MODEL'

   Source for rotational angular frequency :math:`\Orot` values; one
   of

   - :nml:value:`'MODEL'` : Differential rotation, with a spatially varying :math:`\Orot`
     obtained from the stellar model
   - :nml:value:`'UNIFORM'` : Uniform rotation, with a spatially
     constant :math:`\Orot` set by :nml:option:`Omega_rot` and
     :nml:option:`Omega_rot_units`

.. nml:option:: Omega_rot
   :type: real
   :default: 0

   Rotation angular frequency. Used only when :nml:option:`Omega_rot_source`\ =\ :nml:value:`'UNIFORM'`

.. nml:option:: Omega_rot_units
   :type: string
   :default: 'NONE'

   Units of :nml:option:`Omega_rot`; one of

   - :nml:value:`'NONE'` : Dimensionless angular frequency
   - :nml:value:`'HZ'` : Linear frequency in Hz\ [#only-D]_
   - :nml:value:`'UHZ'` : Linear frequency in :math:`\mu`\ Hz\ [#only-D]_
   - :nml:value:`'RAD_PER_SEC'` : Angular frequency in radians per second\ [#only-D]_
   - :nml:value:`'CYC_PER_DAY'` : Linear frequency in cycles per day\ [#only-D]_
   - :nml:value:`'CRITICAL'` : Fraction of the Roche critical rate\ [#only-D]_

   Used only when :nml:option:`Omega_rot_source`\ =\ :nml:value:`'UNIFORM'`

.. nml:option:: tag_list
   :type: string
   :default: ''

   Comma-separated list of :nml:option:`tag <mode.tag>` values to
   match; matches all if left blank

.. rubric:: Footnotes

.. [#only-D] This option is available only for stellar models with :ref:`D capability <model-caps>`

