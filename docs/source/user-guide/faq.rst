.. _faq:

**************************
Frequently Asked Questions
**************************

This is a list of Frequently Asked Questions about GYRE. Suggestions for new entries are always most welcome!

How Do I...
===========

...obtain the GYRE source code?
  See the :ref:`install-download` section.

...compile GYRE?
  See the :ref:`install-compile` section.

.. _faq-multicore:

...run GYRE on multiple cores?
  GYRE can take advantage of multiple
  processors in a shared-memory (multicore) computer through its use
  of :wiki:`OpenMP <OpenMP>`. This functionality should be enabled by
  default, but you can nevertheless force it by setting the :envvar:`OMP`
  environment variable to ``yes`` prior to compilation. Then, set the
  :envvar:`OMP_NUM_THREADS` environment variable to the number of threads
  you want to use.

...pronounce GYRE?
  With a soft 'g' rhyming with 'wire', like :download:`this <faq/gyre-spoken.mp3>`.

...cite GYRE?
  See the :ref:`citing-gyre` section.

...get support for a problem I'm having?
  Post a message to one of the `GYRE discussion forums <gyre-forums_>`__.

...contact the GYRE team?
  Send an email to `Rich Townsend <mailto:townsend@astro.wisc.edu>`__.

...access the documentation for older releases?
  Click on the 'v:...'  dropdown menu in the bottom-right. Note that this
  menu is only available when viewing the
  `documentation <https://gyre.readthedocs.io/en/stable/>`__ on `Read
  the Docs <https://readthedocs.org/>`__; it isn't available in a
  local build within the :file:`{$GYRE_DIR}/docs` directory.
  Also, the legacy documentation is a work-in-progress, and remains
  incomplete.


What Does...
============
...'Failed during deflate narrow : out-of-domain frequency' mean?
  This is an indication that GYRE's root solver wandered out of bounds
  when trying to find a complex root of the discriminant function. Try running
  with a different choice of :nml_n:`diff_scheme` parameter
  (:nml_v:`MAGNUS_GL2` seems to be the most robust), and/or using
  contour method instead (see the :ref:`non-ad-osc` chapter).

Why Does...
===========

...the error 'Illegal Instruction' arise on MacOS when running with large grid sizes?
  This behavior is typically caused by overflow of the OpenMP stack
  (for more info see `here <https://stackoverflow.com/questions/13870564/gfortran-openmp-segmentation-fault-occurs-on-basic-do-loop>`__).
  Try setting the :envvar:`OMP_STACKSIZE` environment variable to 500K or 1M.

