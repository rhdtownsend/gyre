.. _installation:

************
Installation
************

This chapter discusses GYRE installation in detail. If you just want
to get up and running, have a look at the :ref:`quick-start` chapter.

Pre-Requisites
==============

To compile and run GYRE, you'll need the following software
components:

* A modern (2003+) Fortran compiler
* The :netlib:`BLAS <blas>` linear algebra library
* The :netlib:`LAPACK <lapack>` linear algebra library
* The :netlib:`LAPACK95 <lapack95>` Fortran 95
  interface to LAPACK
* The `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`__ data management library
* The `crlibm <https://hal-ens-lyon.archives-ouvertes.fr/ensl-01529804>`__ correctly rounded math library
* The :git:`crmath <rhdtownsend/crmath>` Fortran 2003 interface to crlibm
* An OpenMP-aware version of the :netlib:`ODEPACK <odepack>` differential equation library (optional)

On Linux and MacOS platforms, these components are bundled together in
the MESA Software Development Kit (SDK), which can be downloaded from
the `MESA SDK <mesa-sdk_>`__ homepage. Using this SDK is strongly
recommended.

Building GYRE
=============

.. _install-download:

Download
--------

Download the `GYRE source code <tarball_>`__, and unpack it
from the command line using the :command:`tar` utility:

.. prompt:: bash
   :substitutions:

   tar xf gyre-|version|.tar.gz

Set the :envvar:`GYRE_DIR` environment variable with the path to the
newly created source directory; this can be achieved, e.g., using the
:command:`realpath` command\ [#realpath]_:

.. prompt:: bash
   :substitutions:

   export GYRE_DIR=$(realpath gyre-|version|)

.. _install-compile:

Compile
-------

Compile GYRE using the :command:`make` utility:

.. prompt:: bash

   make -j -C $GYRE_DIR

(the :command:`-j` flags tells :command:`make` to use multiple cores, speeding up the build).

Test
----

To check that GYRE has compiled correctly and gives reasonable
results, you can run the calculation test suite via the command

.. prompt:: bash

   make -C $GYRE_DIR test

The initial output from the tests should look something like this:

.. literalinclude:: installation/make-test.out
   :language: console
   :lines: 1-10

If things go awry, consult the :ref:`troubleshooting`
chapter.

Custom Builds
=============

Custom builds of GYRE can be created by setting certain environment
variables, and/or variables in the file
:file:`{$GYRE_DIR}/src/build/Makefile`, to the value ``yes``. The
following variables are currently supported:

DEBUG
  Enable debugging mode (default ``no``)

OMP
  Enable OpenMP parallelization (default ``yes``)

MPI
  Enable MPI parallelization (default ``no``)

DOUBLE_PRECISION
  Use double precision floating point arithmetic (default ``yes``)

CRMATH
  Use correctly rounded math functions (default ``yes``)

IEEE
  Use Fortran IEEE floating point features (default ``no``)

FPE
  Enable floating point exception checks (default ``yes``)

HDF5
  Include HDF5 support (default ``yes``)

EXPERIMENTAL
  Enable experimental features (default ``no``)

If a variable is not set, then its default value is assumed.

Git Access
==========

Sometimes, you'll want to try out new features in GYRE that haven't
yet made it into a formal release. In such cases, you can check out
GYRE directly from the :git:`rhdtownsend/gyre` git repository on
:git:`GitHub <>`:

.. prompt:: bash

   git clone https://github.com/rhdtownsend/gyre.git

However, a word of caution: GYRE is under constant development, and
features in the main (``master``) branch can change without warning.

.. rubric:: footnote

.. [#realpath] The :command:`realpath` command is included in the GNU
               `CoreUtils <https://www.gnu.org/software/coreutils/>`__
               package. Mac OS users can install CoreUtils using
               `MacPorts <https://www.macports.org/>`__ or `Homebrew
               <https://brew.sh/>`__.
