.. _model-group:

.. nml:group:: model

Model Namelist Group
====================

The :nml:group:`model` namelist group controls how the stellar model is
set up. The following options are available:

.. nml:option:: model_type
   :type: string
   :default: ''

   Type of stellar model; one of

   - :nml:value:`'HOM'` : :ref:`Homogeneous compressible model <hom-models>`
   - :nml:value:`'POLY'` : :ref:`Polytropic model <poly-models>` read from external file
   - :nml:value:`'ANAPOLY_0'` : :ref:`Analytic polytropic model <anapoly-models>` with :math:`\npoly=0`
   - :nml:value:`'ANAPOLY_1'` : :ref:`Analytic polytropic model <anapoly-models>` with :math:`\npoly=1`
   - :nml:value:`'ANAPOLY_5'` : :ref:`Analytic polytropic model <anapoly-models>` with :math:`\npoly=5`
   - :nml:value:`'ANAPOLY_5_1'` : :ref:`Analytic polytropic model <anapoly-models>` with :math:`\npoly=5,1`
   - :nml:value:`'EVOL'` : :ref:`Evolutionary model <evol-models>` read from external file

.. nml:option:: file
   :type: string
   :default: ''

   Name of file. Used only when :nml:option:`model_type`\ = \
   :nml:valuelist:`'POLY' 'EVOL'`

.. nml:option:: file_format
   :type: string
   :default: ''

   Format of file; one of

   - :nml:value:`'AMDL'` : AMDL-format binary file
   - :nml:value:`'B3'` : B3-format HDF5 file
   - :nml:value:`'FAMDL'` : FAMDL-format text file
   - :nml:value:`'FGONG'` : FGONG-format text file
   - :nml:value:`'GSM'` : :ref:`GSM-format <gsm-file-format>` HDF5 file
   - :nml:value:`'LOSC'` : LOSC-format text file
   - :nml:value:`'MESA'` : :ref:`MESA/GYRE-format <mesa-file-format>` text file
   - :nml:value:`'OSC'` : OSC-format text file
   - :nml:value:`'WDEC'` : WDEC-format text file

   Used only when :nml:option:`model_type`\ = \ :nml:value:`'EVOL'`

.. nml:option:: data_format
   :type: string
   :default: ''

   Fortran format specifier for data read from OSC-, FGONG- and
   FAMDL-format files. If left blank, format is auto-selected

.. nml:option:: deriv_type
   :type: string
   :default: 'MONO'

   Cubic interpolation derivatives type; one of

   - :nml:value:`'SPLINE'` : Spline (non-local) derivatives
   - :nml:value:`'FINDIFF'` : Finite-difference derivatives
   - :nml:value:`'MONO'` : Monotonized derivatives

   Used only when :nml:option:`model_type` = :nml:value:`'EVOL'` and
   :nml:option:`interp_type` = :nml:value:`'CUBIC'`

.. nml:option:: Gamma_1
   :type: real
   :default: 5/3

   First adiabatic exponent. Used only when :nml:option:`model_type`
   = :nml:valuelist:`'HOM' 'ANAPOLY_0' 'ANAPOLY_1' 'ANAPOLY_5'
   'ANAPOLY_5_1'`

.. nml:option:: theta_s
   :type: real
   :default: 0

   Surface value of polytropic dependent variable. Used only when
   :nml:option:`model_type` = :nml:valuelist:`'ANAPOLY_0'
   'ANAPOLY_1' 'ANAPOLY_5'`

.. nml:option:: x_match
   :type: real
   :default: 0.5

   Radial coordinate of match point between inner and outer regions.
   Used only when :nml:option:`model_type` =
   :nml:value:`'ANAPOLY_5_1'`

.. nml:option:: grid_type
   :type: string
   :default: 'UNI'

   Model grid type; one of

   - :nml:value:`'UNI'` : Uniform spacing
   - :nml:value:`'GEO'` : Geometric spacing
   - :nml:value:`'LOG'` : Logarithmic spacing

   Used only when :nml:option:`model_type` = :nml:valuelist:`'HOM' 'ANAPOLY_0' 'ANAPOLY_1' 'ANAPOLY_5' 'ANAPOLY_5_1'`

.. nml:option:: n
   :type: integer
   :default: 10

   Number of points in model grid. Used only when
   :nml:option:`model_type` = :nml:valuelist:`'HOM' 'ANAPOLY_0' 'ANAPOLY_1' 'ANAPOLY_5' 'ANAPOLY_5_1'`

.. nml:option:: s
   :type: real
   :default: 1

   Skewness coefficient for model grid. Used only when
   :nml:option:`model_type` = :nml:valuelist:`'HOM' 'ANAPOLY_0'
   'ANAPOLY_1' 'ANAPOLY_5' 'ANAPOLY_5_1'` and :nml:option:`grid_type`
   = :nml:valuelist:`'GEO LOG'`

.. nml:option:: x_i
   :type: real
   :default: 0

   Inner boundary coordinate of model grid. Used only when
   :nml:option:`model_type` = :nml:valuelist:`'HOM' 'ANAPOLY_0'
   'ANAPOLY_1' 'ANAPOLY_5' 'ANAPOLY_5_1'`

.. nml:option:: x_o
   :type: real
   :default: 1

   Outer boundary coordinate of model grid. Used only when
   :nml:option:`model_type` = :nml:valuelist:`'HOM' 'ANAPOLY_0'
   'ANAPOLY_1' 'ANAPOLY_5' 'ANAPOLY_5_1'`

.. nml:option:: dx_snap
   :type: real
   :default: 0

   Threshold for snapping model points together; if a pair of points
   are separated by less than :nml:option:`dx_snap`, they are snapped
   together. Used only when :nml:option:`model_type` =
   :nml:value:`'EVOL'`

.. nml:option:: add_center
   :type: logical
   :default: .TRUE.

   Flag to add a center point to the model; if a point does not
   already exist at the origin, then one is added. Used only when
   :nml:option:`model_type` = :nml:valuelist:`'EVOL' 'POLY'`

.. nml:option:: interp_type
   :type: string
   :default: 'CUBIC'

   Interpolation type; one of

   - :nml:value:`'CUBIC'`  : Piecewise cubic
   - :nml:value:`'LINEAR'` : Piecewise linear

   Used only when :nml:option:`model_type` = :nml:valuelist:`'EVOL' 'POLY'`

.. nml:option:: repair_As
   :type: logical
   :default: .FALSE.

   Flag to repair inaccuracies in the dimensionless Brunt-Väisälä
   frequency at density discontinuities

.. nml:option:: constrain_derivs
   :type: logical
   :default: .TRUE.

   Flag to constrain first derivatives of :math:`V_2`, :math:`U` and
   :math:`c_1` structure coefficients, in accordance with equations
   (20) and (21) of :ads_citet:`takata:2006a` and the hydrostatic
   equilibrium equation. Used only when :nml:option:`model_type` =
   :nml:value:`'EVOL' 'POLY'`

.. nml:option:: use_nabla_rad
   :type: logical
   :default: .FALSE.

   Obtain the radiative luminosity :math:`\Lrad` from :math:`\nabrad`. Used
   only when :nml:option:`model_type` = :nml:value:`'EVOL'` and
   :nml:option:`file_format` = :nml:value:`'OSC'`
