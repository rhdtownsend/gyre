.. _build-poly-ex-comp:

Example Walkthrough: Composite Polytrope
========================================

As the second example of :program:`build_poly` in action, let's build
a two-region composite polytrope. The polytropic index is :math:`n=3`
in the inner region, and :math:`n=1.5` in the outer region. At the
boundary between the regions, located at radial coordinate
:math:`z=1.4`, the logarithmic density jump is :math:`\Delta = -0.5`.

Assembling a Namelist File
--------------------------

Using a text editor, create the file :file:`build_poly.composite.in` with
the following content cut-and-pasted in:

.. literalinclude:: build_poly.composite.in

Again, detailed information on the namelist groups expected in
:program:`build_poly` input files can be found in the
:ref:`build-poly-input` section. Here, let's briefly narrate the
parameters appearing in the file above:

* In the :nml_g:`poly` namelist group, the :nml_n:`n_r` parameter sets
  the number of regions; the :nml_n:`n_poly` parameter sets the
  polytropic indices in the two regions; the :nml_n:`z_b` sets the
  radial coordinate of the boundary between the regions; and the
  :nml_n:`Delta_b` sets the density jump at this boundary.
* In the :nml_g:`num` namelist group, the :nml_n:`dz` parameter sets
  the radial spacing of points, and the :nml_n:`toler` parameter sets
  the tolerance of the numerical integrator.
* In the :nml_g:`output` namelist group, the :nml_n:`file` parameter
  sets the name of the output file.

Running build_poly
------------------

As before, to run :program:`build_poly` use the command

.. code-block:: console
			 
   $ $GYRE_DIR/bin/build_poly build_poly.composite.in

There is no screen output produced during the run, but at the end the
:file:`poly.composite.h5` will be written to disk. This file, which is in
:ref:`POLY format<poly-file-format>`, can be used as the input stellar
model in a GYRE calculation; but it can also be explored in Python
(see :numref:`fig-poly-comp`) using the ``read_model`` function from
:git:`PyGYRE <rhdtownsend/pygyre>`.

.. _fig-poly-comp:

.. figure:: fig_poly_composite.svg
   :alt: Plot showing the structure of the simple polytrope model
   :align: center

   Plot of the Lane-Emden solution variable :math:`\theta`, density
   :math:`\rho`, pressure :math:`P` and interior mass :math:`M_{r}` as a
   function of radial coordinate, for the composite polytrope. Note
   the density discontinuity, and the associated discontinuities in
   the gradients of the pressure and interior mass. (:download:`Source
   <fig_poly_composite.py>`)
