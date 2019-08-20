Mode Parameters
===============

The :nml_g:`mode` namelist group defines mode parameters; the input
file can contain one or more. Allowable fields are:

:nml_o:`l` (default :nml_l:`0`)
  Harmonic degree :math:`\ell`
  
:nml_o:`m` (default :nml_l:`0`)
  Azimuthal order :math:`m`

:nml_o:`tag`
  Tag for controlling selection of other parameters

:nml_o:`n_pg_min` (default :nml_l:`-HUGE`)
  Filter for minimum radial order

:nml_o:`n_pg_max` (default :nml_l:`+HUGE`)
  filter for maximum radial order

:nml_o:`rossby` (default :nml_l:`.FALSE.`)
  Flag to use Rossby-mode angular eigenvalues/eigenfunctions
