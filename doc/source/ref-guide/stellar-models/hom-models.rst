.. _hom-models:

Homogeneous Models
==================

Setting the :nml_n:`model_type` parameter of the :nml_g:`model`
namelist group to :nml_v:`'HOM'` tells the frontend to create a
homogeneous (uniform density) stelllar model, equivalent to a
polytrope with index :math:`n=0`. Because the structure of these model
can be computed analytically, there is no need to read from an
external file.

The :nml_n:`Gamma_1` parameter of the :nml_g:`model` namelist group
controls the first adiabatic index of the model, while the :nml_n:`n`,
:nml_n:`s` and :nml_n:`grid_type` parameters control the model
grid. See the :ref:`model-params` section for further details.
