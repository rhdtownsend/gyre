Output Parameters
=================

The :nml_g:`ad_output` and :nml_g:`nad_output` namelist groups
determine the output produced at the end of a run, from the adiabatic
and non-adiabatic calculation stages, respectively; the input file
should contain exactly one of each. Allowable fields are:

:nml_o:`summary_file` (default :nml_l:`''`)
  Name of summary file

:nml_o:`summary_file_format` (default :nml_l:`'HDF'`)
  Format of summary file; one of

  :nml_l:`'HDF'` : HDF5 file
  :nml_l:`'TXT'` : Text file

:nml_o:`summary_item_list` (default :nml_l:`''`)
  Comma-separated list of output items to write to summary file; see `Output
  Files <Output Files (5.1)>`__ for possible choices

:nml_o:`mode_template` (default :nml_l:`''`)
  Name template of mode files. Names are generated using the following pattern
  substitutions:

  - :nml_l:`'%J'` : Unique mode index :math:`j`, formatted in fixed-width field
  - :nml_l:`%j` : Same as ``%J``, but formatted in variable-width field
  - :nml_l:`%L` : Harmonic degree :math:`\ell`, formatted in fixed-width field
  - :nml_l:`%l` : Same as ``%L``, but formatted in variable-width field
  - :nml_l:`%N` : Radial order :math:`n_{\rm pg}`, formatted in fixed-width field
  - :nml_l:`%n` : Same as ``%N``, but formatted in variable-width field

:nml_o:`mode_file_format` (default :nml_l:`'HDF'`)
  Format of mode files; one of

  - :nml_o:`'HDF'` : HDF5 file
  -  :nml_l:`'TXT'` : text file

:nml_o:`mode_item_list` (default :nml_l:`''`)
  Comma-separated list of output items to write to mode files; see `Output
  Files <Output Files (5.1)>`__ for possible choices

:nml_o:`freq_units` (default :nml_l:`NONE`)
  Units of :nml_l:`freq` output item; one of:

  - :nml_l:`'NONE'` : Dimensionless angular frequency
  - :nml_l:`'HZ'` : linear frequency in Hz [#only_evol]_
  - :nml_l:`'UHZ'` : linear frequency in μHz [#only_evol]_
  - :nml_l:`'RAD_PER_SEC'` : angular frequency in radians per second [#only_evol]_
  - :nml_l:`'CYC_PER_DAY'` : linear frequency in cycles per day [#only_evol]_
  - :nml_l:`'ACOUSTIC_DELTA'` : Fraction of the asymptotic acoustic large frequency separation :math:`\Delta \nu`
  - :nml_l:`'GRAVITY_DELTA'` : Fraction of the asymptotic inverse gravity period separation :math:`(\Delta P)^{-1}`
  - :nml_l:`'UPPER_DELTA'` : Greater of :math:`\Delta \nu` and :math:`(\Delta P)^{-1}`
  - :nml_l:`'LOWER_DELTA'` : Lesser of :math:`\Delta \nu` and :math:`(\Delta P)^{-1}`
  - :nml_l:`'ACOUSTIC_CUTOFF'` : fraction of the acoustic cutoff frequency [#only_evol]_
  - :nml_l:`'GRAVITY_CUTOFF'` : fraction of the gravity cutoff frequency [#only_evol]_
  - :nml_l:`'ROSSBY_I'` : fraction of Rossby frequency at inner boundary
  - :nml_l:`'ROSSBY_O'` : fraction of Rossby frequency at outer boundary

:nml_o:`freq_frame`` (default :nml_l:`INERTIAL`)
  Frame of :nml_l:`freq` output item; one of:

   - :nml_l:`'INERTIAL'` : Inertial frame
   - :nml_l:`'COROT_I'` : Co-rotating frame at inner boundary
   - :nml_l:`'COROT_O'` : Co-rotating frame at outer boundary

:nml_o:`label` (default :nml_l:`''`)
  Textual label to add to all output files

:nml_o:`prune_modes` (default :nml_l:`.FALSE.`)
  Flag to discard eigenfunction data after (possibly) writing it to
  disk; used to conserve memory

.. rubric:: Footnotes

.. [#only_evol] This option is only available when :nml_o:`model_type` is :nml_l:`'EVOL'`
