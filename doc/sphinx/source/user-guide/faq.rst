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
  environment variable to `yes` prior to compilation. Then, set the
  :envvar:`OMP_NUM_THREADS` environment variable to the number of threads
  you want to use.

.. _faq-cluster:

...run GYRE on a cluster?
  GYRE can take advantage of multiple nodes in a computer cluster
  through its use of :wiki:`MPI <Message Passing Interface>`. To enable this functionality,
  set the :envvar:`MPI` environment variable to `yes` prior to
  compilation. Note that you'll need to use a version of the
  `MESA Software Development Kit <mesa-sdk_>`__ built with MPI
  support for your cluster's specific networking technology (e.g.,
  :wiki:`Infiniband <Infiniband>`); contact the GYRE team for more details.

...pronounce GYRE?
  With a soft 'g' rhyming with 'wire', like :download:`this <faq/gyre-spoken.mp3>`.

...cite GYRE?
  See the :ref:`citing-gyre` section.

...get support for a problem I'm having?
  Post a message to one of the `GYRE discussion forums <gyre-forums_>`__.

...access the documentation for older releases?
  Click on the 'v:...'  dropdown menu at the bottom of the left-hand
  panel. Note that the legacy documentation is a work-in-progress,
  and remains incomplete.


What Does...
============
...'Failed during deflate narrow : out-of-domain frequency' mean?
  This is an indication that GYRE's root solver wandered out of bounds
  when trying to find a complex root of the discriminant function. Try running
  with a different choice of :nml_n:`diff_scheme` parameter
  (:nml_v:`MAGNUS_GL2` seems to be the most robust), and/or using
  :program:`gyre_contour` instead.

Why Does...
===========

...the error 'Illegal Instruction' arise on MacOS when running with large grid sizes?
  This behavior is typically caused by overflow of the OpenMP stack
  (for more info see `here <http://stackoverflow.com/questions/13870564/gfortran-openmp-segmentation-fault-occurs-on-basic-do-loop>`__).
  Try setting the :envvar:`OMP_STACKSIZE` environment variable to 500K or 1M.

