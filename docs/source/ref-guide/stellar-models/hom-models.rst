.. _hom-models:

Homogeneous Models
==================

.. nml:group:: model
   :no-target:

Setting the :nml:option:`model_type` option of the :nml:group:`model`
namelist group to :nml:value:`'HOM'` tells the frontend to create a
homogeneous (uniform density) stellar model, equivalent to a
polytrope with index :math:`n=0`. Because the structure of these model
can be computed analytically, there is no need to read from an
external file.

The :nml:option:`Gamma_1` option controls the first adiabatic index of
the model, while the :nml:option:`n`, :nml:option:`s` and
:nml:option:`grid_type` options control the model grid. See the
:ref:`model-group` section for further details.
