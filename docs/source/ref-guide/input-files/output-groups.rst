.. _output-groups:

.. nml:group:: ad_output
.. nml:group:: nad_output
.. nml:group:: tide_output
.. nml:group:: output

Output Namelist Groups
======================

The :nml:group:`ad_output`, :nml:group:`nad_output` and
:nml:group:`tide_output` namelist groups determine the output
produced at the end of a run (the first two for the adiabatic and
non-adiabatic calculation stages of :program:`gyre`; the third for
:program:`gyre_tides`). The following options are available:

.. nml:option:: summary_file
   :type: string
   :default: ''

   Name of summary file; if left blank, no summary file is written

.. nml:option:: summary_file_format
   :type: string
   :default: 'HDF'

   Format of summary file; one of

   - :nml:value:`'HDF'` : HDF5 file
   - :nml:value:`'TXT'` : Text file

.. nml:option:: summary_item_list
   :type: string
   :default: 'l,n_pg,omega,freq'

   Comma-separated list of fields to write to summary file; see
   the :ref:`summary-files` section for possible choices

.. nml:option:: summary_filter_list
   :type: string
   :default: ''

   Comma-separated list of filter criteria for summary files; see the
   :ref:`output-filters` section for possible choices

.. nml:option:: detail_template
   :type: string
   :default: ''

   Name template of detail files. If left blank, no detail files are
   written. Names are generated from the template by applying the
   following pattern substitutions:

   - :nml:value:`'%ID'` : Unique mode index, formatted in fixed-width field
   - :nml:value:`'%id'` : Same as :nml:value:`'%ID'`, but formatted in variable-width field
   - :nml:value:`'%L'` : Harmonic degree :math:`\ell`, formatted in fixed-width field
   - :nml:value:`'%l'` : Same as :nml:value:`'%L'`, but formatted in variable-width field
   - :nml:value:`'%M'` : Azimuthal order :math:`m`, formatted in fixed-width field
   - :nml:value:`'%m'` : Same as :nml:value:`'%M'`, but formatted in variable-width field
   - :nml:value:`'%N'` : Radial order :math:`n_{\rm pg}`, formatted in fixed-width field
   - :nml:value:`'%n'` : Same as :nml:value:`'%N'`, but formatted in variable-width field
   - :nml:value:`'%P'` : Acoustic wave winding number :math:`n_{\rm p}`, formatted in fixed-width field
   - :nml:value:`'%p'` : Same as :nml:value:`'%P'`, but formatted in variable-width field
   - :nml:value:`'%G'` : Gravity wave winding number :math:`n_{\rm g}`, formatted in fixed-width field
   - :nml:value:`'%g'` : Same as :nml:value:`'%G'`, but formatted in variable-width field

.. nml:option:: detail_file_format
   :type: string
   :default: 'HDF'

   Format of detail files; one of

   - :nml:value:`'HDF'` : HDF5 file
   - :nml:value:`'TXT'` : Text file

.. nml:option:: detail_item_list
   :type: string
   :default: 'l,n_pg,omega,freq,x,xi_r,xi_h'

   Comma-separated list of fields to write to detail files; see the
   :ref:`detail-files` section for possible choices

.. nml:option:: detail_filter_list
   :type: string
   :default: ''

   Comma-separated list of filter criteria for detail files; see the
   :ref:`output-filters` section for possible choices

.. nml:option:: freq_units
   :type: string
   :default: 'NONE'

   Units of frequency-like output fields; one of

   - :nml:value:`'NONE'` : Dimensionless angular frequency
   - :nml:value:`'HZ'` : Linear frequency in Hz\ [#only-D]_
   - :nml:value:`'UHZ'` : Linear frequency in :math:`\mu`\ Hz\ [#only-D]_
   - :nml:value:`'RAD_PER_SEC'` : Angular frequency in radians per second\ [#only-D]_
   - :nml:value:`'CYC_PER_DAY'` : Linear frequency in cycles per day\ [#only-D]_
   - :nml:value:`'ACOUSTIC_DELTA'` : Fraction of the asymptotic acoustic large frequency separation :math:`\Delta \nu`
   - :nml:value:`'GRAVITY_DELTA'` : Fraction of the asymptotic inverse gravity period separation :math:`(\Delta P)^{-1}`
   - :nml:value:`'UPPER_DELTA'` : Greater of :math:`\Delta \nu` and :math:`(\Delta P)^{-1}`
   - :nml:value:`'LOWER_DELTA'` : Lesser of :math:`\Delta \nu` and :math:`(\Delta P)^{-1}`
   - :nml:value:`'ACOUSTIC_CUTOFF'` : Fraction of the acoustic cutoff frequency\ [#only-D]_
   - :nml:value:`'GRAVITY_CUTOFF'` : Fraction of the gravity cutoff frequency\ [#only-D]_
   - :nml:value:`'ROSSBY_I'` : Fraction of Rossby frequency at inner boundary
   - :nml:value:`'ROSSBY_O'` : Fraction of Rossby frequency at outer boundary

.. nml:option:: freq_frame
   :type: string
   :default: 'INERTIAL'`

   Frame of frequency-like output fields; one of

   - :nml:value:`'INERTIAL'` : Inertial frame
   - :nml:value:`'COROT_I'` : Co-rotating frame at inner boundary
   - :nml:value:`'COROT_O'` : Co-rotating frame at outer boundary

.. nml:option:: label
   :type: string
   :default: ''

   Textual label to add to all output files

.. rubric:: Footnotes

.. [#only-D] This option is available only for stellar models with :ref:`D capability <model-caps>`
