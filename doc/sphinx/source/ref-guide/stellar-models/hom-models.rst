.. _hom-models:

Homogeneous Models
==================

Homogeneous models represent uniform-density stars in hydrostatic
equilibrium; they are equivalent to polytropes with index
:math:`\npoly=0`. Because their structure can be calculated
analytically, GYRE creates them on-the-fly without the need to read
from an external file.

To use a homogeneous model with GYRE, set the :nml_n:`model_type`
parameter in the :nml_g:`model` namelist group to :nml_v:`'HOM'`, and
set the :nml_n:`Gamma_1`, :nml_n:`n`, :nml_n:`s` and
:nml_n:`grid_type` parameters to suitable values (see the
:ref:`model-params` section for more details).
