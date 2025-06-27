.. _build-poly-ex-simple:

Example Walkthrough: Simple Polytrope
=====================================

As the first example of :program:`build_poly` in action, let's build a
simple (i.e., single-region) :math:`n=3` polytrope, that for instance
describes the structure of a radiation-pressure dominated, fully
convective star.

Assembling a Namelist File
--------------------------

First, let's assemble a namelist file containing the various
parameters which control a :program:`build_poly` run.  Using a text
editor, create the file :file:`build_poly.simple.in` with the
following content cut-and-pasted in:

.. literalinclude:: build_poly.simple.in

Detailed information on the namelist groups expected in
:program:`build_poly` input files can be found in the
:ref:`build-poly-input` section. Here, let's briefly narrate the
parameters appearing in the file above:

* In the :nml_g:`poly` namelist group, the :nml_n:`n_poly` parameter
  sets the polytropic index.
* In the :nml_g:`num` namelist group, the :nml_n:`dz` parameter sets
  the radial spacing of points, and the :nml_n:`toler` parameter sets
  the tolerance of the numerical integrator.
* In the :nml_g:`output` namelist group, the :nml_n:`file` parameter
  sets the name of the output file.

Running build_poly
------------------

To run :program:`build_poly`, use the command

.. code-block:: console
			 
   $ $GYRE_DIR/bin/build_poly build_poly.simple.in

There is no screen output produced during the run, but at the end the
:file:`poly.simple.h5` will be written to disk. This file, which is in
:ref:`POLY format<poly-file-format>`, can be used as the input stellar
model in a GYRE calculation; but it can also be explored in Python
(see :numref:`fig-poly-simple`) using the ``read_model`` function from
:git:`PyGYRE <rhdtownsend/pygyre>`.

.. _fig-poly-simple:

.. figure:: fig_poly_simple.svg
   :alt: Plot showing the structure of the simple polytrope model
   :align: center

   Plot of the Lane-Emden solution variable :math:`\theta`, density
   :math:`\rho`, pressure :math:`P` and interior mass :math:`M_{r}` as a
   function of radial coordinate, for the simple
   polytrope. (:download:`Source <fig_poly_simple.py>`)
