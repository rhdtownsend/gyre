.. _hom-models:

Homogeneous Models
==================

Homogeneous models represent a uniform-density star in hydrostatic
equilibrium; they are equivalent to a polytrope with index
:math:`\npoly=0`. Because their structure can be calculated
analytically, GYRE constructs them on-the-fly without the need to read
from a file. The model construction is controlled by the
:nml_n:`Gamma_1`, :nml_n:`n`, :nml_n:`s` and :nml_n:`grid_type`
parameters of the :nml_g:`model` namelist group (see the
:ref:`model-params` section for more information).

Although homogeneous models are a poor approximation to real stars,
they are very useful for testing purposes because analytic expressions
exist for their eigenfrequencies and eigenfunctions (see, e.g.,
:ads_citealp:`pekeris:1938`).
