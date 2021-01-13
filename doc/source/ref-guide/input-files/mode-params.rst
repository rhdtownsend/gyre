.. _mode-params:

Mode Parameters
===============

The :nml_g:`mode` namelist group defines mode parameters; the input
file can contain one or more. Allowable parameters are:

:nml_n:`l` (default :nml_v:`0`)
  Harmonic degree :math:`\ell`
  
:nml_n:`m` (default :nml_v:`0`)
  Azimuthal order :math:`m`

:nml_n:`tag`
  Tag for controlling selection of other parameters

:nml_n:`n_pg_min` (default :nml_v:`-HUGE`)
  Filter for minimum radial order (applies only to adiabatic calculations)

:nml_n:`n_pg_max` (default :nml_v:`+HUGE`)
  Filter for maximum radial order (applies only to adiabatic calculations)

:nml_n:`static` (default :nml_v:`.FALSE.`)
  Flag to solve for the static (:math:`\omega \rightarrow 0`) limit
