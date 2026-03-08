.. _walkthrough:

*******************
Example Walkthrough
*******************

This chapter provides a walkthrough of a example GYRE project, to
illustrate the typical steps involved. For this example, we'll use
:program:`gyre` (the :ref:`frontend <frontends>` focused on stellar
oscillations) to find eigenfrequencies and eigenfunctions of dipole
and quadrupole gravity modes for a MESA model of slowly pulsating B
(SPB) star.

.. _walkthrough-work:

Making a Place to Work
======================

When starting a new project, it's a good idea to create a dedicated
work directory to contain the various input and output files that :program:`gyre`
operates on. These commands will make a new directory beneath your
home directory with the name :file:`work`, and then set this directory
as the current working directory:

.. code-block:: console

   $ mkdir ~/work
   $ cd ~/work

Grabbing a Stellar Model
========================

The next step is to grab the stellar model. There are a number of
example models provided in the :file:`${GYRE_DIR}/models` directory;
the following commands will copy a MESA model for a :math:`5\,\Msun`
SPB star into your work directory:

.. code-block:: console

   $ cp $GYRE_DIR/models/mesa/spb/spb.mesa .

Assembling a Namelist File
==========================

Now comes the fun part: assembling an input file containing the
various settings that control a :program:`gyre` run. Using a text
editor, create the file :file:`gyre.in` in your work directory with
the following content copy-and-pasted in (use the copy button that
appears when you hover over the top-right corner of the box):

.. literalinclude:: example-walkthrough/gyre.in

This file is in namelist format, containing multiple namelist
groups. Detailed information on the groups can be found in the
:ref:`namelist-input-files` chapter; for now, let's just focus on some
of the more-important aspects of the file above:

* the :nml:group:`constants` namelist group is empty, telling
  :program:`gyre` to use default values for fundamental constants;
* the :nml:group:`model` namelist group tells :program:`gyre` to read
  an evolutionary model, in :ref:`MESA format <mesa-file-format>`,
  from the file :file:`spb.mesa`;
* the two :nml:group:`mode` namelist groups tells :program:`gyre` to
  search first for dipole (:math:`\ell=1`) and then quadrupole
  (:math:`\ell=2`) modes;
* the :nml:group:`osc` namelist group tells :program:`gyre` to assume,
  when setting up the outer boundary conditions in the oscillation
  equations, that the density vanishes at the stellar surface;
* the :nml:group:`scan` namelist group tells :program:`gyre` to scan a
  region of dimensionless angular frequency space typically occupied
  by gravity modes;
* the :nml:group:`grid` namelist group tells :program:`gyre` how to
  refine the model spatial grid;
* the :nml:group:`ad_output` namelist group tells :program:`gyre` what
  adiabatic data to write to which output files; summary data to the
  file :file:`summary.h5`, and individual mode data to files having
  the prefix ``mode.``;
* the :nml:group:`nad_output` namelist group is empty, telling
  :program:`gyre` not to write any non-adiabatic data.

Running gyre
============

With the hard work done, it's now trivial to run :program:`gyre`:

.. code-block:: console

   $ $GYRE_DIR/bin/gyre gyre.in

As the frontend runs (on multiple cores, if you have a multi-core machine;
see the :ref:`FAQ <faq-multicore>` for more details), it will print lots of data
to the screen. Let's break down this output, chunk by chunk.

First, :program:`gyre` prints out its version number, tells us (in
OpenMP threads) how many cores it is running on, and indicates which
file it is reading settings from (here, :file:`gyre.in`):

.. literalinclude:: example-walkthrough/gyre.out
   :language: console
   :end-before: Model Init

Next, :program:`gyre` loads the stellar model from the file
:file:`spb.mesa`. This model comprises 1814 points and extends from
the surface all the way to the center (which is why :program:`gyre` decides not
to add a central point).

.. literalinclude:: example-walkthrough/gyre.out
   :language: console
   :start-at: Model Init
   :end-before: Mode Search

:program:`gyre` then prepares to search for modes with harmonic degree
:math:`\ell=1` and azimuthal order :math:`m=0` (not specified in
:file:`gyre.in`, but assumed by default), by building a frequency grid
and a spatial grid:

.. literalinclude:: example-walkthrough/gyre.out
   :language: console
   :start-at: Mode Search
   :end-before: Starting search

(The concepts of spatial and frequency grids are explored in greater
detail in the :ref:`numerical` and :ref:`understanding-grids`
chapters). Next, :program:`gyre` attempts to bracket roots of the discriminant
function (again, see the :ref:`numerical` chapter) by
scanning for changes in its sign:

.. literalinclude:: example-walkthrough/gyre.out
   :language: console
   :start-at: Starting search
   :end-before: Finding roots

Finally, for each bracket found :program:`gyre` uses a root solver to
converge to the eigenfrequency, and then processes the mode to
identify it and perform other calculations. Each row of output here
corresponds to a mode that :program:`gyre` has successfully found:

.. literalinclude:: example-walkthrough/gyre.out
   :language: console
   :start-at: Finding roots
   :end-before: Mode Search

The columns appearing are as follows:

.. ofile:filetype:: summary
   :no-target:

:ofile:field:`id`
   Unique mode index

:ofile:field:`l`
   Harmonic degree :math:`\ell`

:ofile:field:`m`
   Azimuthal order :math:`m`

:ofile:field:`n_pg`
   Radial order :math:`\numpg` within the
   Eckart-Scuflaire-Osaki-Takata scheme (see
   :ads_citealp:`takata:2006b`)

:ofile:field:`n_p`
   Acoustic-wave winding number :math:`\nump`

:ofile:field:`n_g`
   Gravity-wave winding number :math:`\numg`

``Re(omega)``
  real part of dimensionless eigenfrequency :math:`\omega`

``Im(omega)``
  imaginary part of dimensionless eigenfrequency :math:`\omega` (zero
  here because we've performed an adiabatic calculation)

:ofile:field:`chi`
   Root-finding convergence parameter :math:`\chi`

:ofile:field:`n_iter`
   Root-finding number of iterations

These values are printed to screen primarily to give an idea of
:program:`gyre`'s progress. Some things to watch out for:

* The convergence parameter :ofile:field:`chi`, defined as the ratio
  of discriminant values before and after the root finding, should
  small (on the order of 1E-9 to 1E-15). If it is significantly larger
  than this, the mode may not be properly converged; and if it is
  significantly smaller than this, there may be numerical issues with
  the discretization scheme.

* The number of iterations :ofile:field:`n_iter` should be moderate;
  values above 20 or so indicate that :program:`gyre` is having
  problems converging.

* The mode radial order :ofile:field:`n_pg` should be
  monotonic-increasing. Departures from this behavior can happen for a
  number of reasons, that are discussed in the :ref:`troubleshooting`
  chapter.

After processing the dipole modes, :program:`gyre` repeats the search
steps for the quadrupole modes. Once the overall run is complete, a
number of output files are written:

* A summary file with the name :file:`summary.h5`

* For each mode found, a detail file with the name
  :file:`detail.l{L}.n{N}.h5`, where :file:`{L}` and :file:`{N}` are
  the harmonic degree and radial order of the mode, respectively.

The :ref:`interpreting-output` chapter discusses how to read and analyze
these files.
