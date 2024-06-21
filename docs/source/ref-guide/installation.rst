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
* A `Python <https://www.python.org/>`__ interpreter (3.5 or more recent)

On Linux and MacOS platforms, all of these components (save for
Python 3) are bundled together in the MESA Software Development Kit
(SDK), which can be downloaded from the `MESA SDK <mesa-sdk_>`__
homepage. Using this SDK is strongly recommended.

Building GYRE
=============

.. _install-download:

Download
--------

Download the `GYRE source code <tarball_url_>`__, and unpack it
from the command line using the :command:`tar` utility:

.. prompt:: bash
   :substitutions:

   tar xf |tarball|

Set the :envvar:`GYRE_DIR` environment variable with the path to the
newly created source directory; this can be achieved, e.g., using the
:command:`realpath` command\ [#realpath]_:

.. prompt:: bash
   :substitutions:

   export GYRE_DIR=$(realpath |dist_dir|)

.. _install-compile:

Compile
-------

Compile GYRE using the :command:`make` utility:

.. prompt:: bash

   make -j -C $GYRE_DIR install

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
variables to the value ``yes``. The following variables are currently
supported:

TOOLS
  Build support tools such as :program:`build_poly` and :program:`eval_lambda` (default ``yes``)

IFACES
  Build additional libraries for interfacing with other codes (default ``no``)

FORUM
  Build the ForUM library internally (default ``yes``)

DEBUG
  Enable debugging mode (default ``no``)

SHARED
  Build shared libraries in addition to static ones (default ``yes``)

OMP
  Enable OpenMP parallelization (default ``yes``)

FPE
  Enable floating point exception checks (default ``yes``)

CRMATH
  Enable correctly rounded math functions (default ``yes``)

PORTABLE
  Enable portable math (default ``yes``)

IEEE
  Use Fortran IEEE floating point features (default ``yes``)

EXPERIMENTAL
  Enable experimental features (default ``no``)

If an environment variable is not set, then its default value is assumed.

Git Access
==========

Sometimes, you'll want to try out new features in GYRE that haven't
yet made it into a formal release. In such cases, you can check out
GYRE directly from the :git:`rhdtownsend/gyre <rhdtownsend/gyre>` git
repository on :git:`GitHub <>`:

.. prompt:: bash

   git clone --recurse-submodules https://github.com/rhdtownsend/gyre.git

However, a word of caution: GYRE is under constant development, and
features in the main (``master``) branch can change without warning.

.. rubric:: footnote

.. [#realpath] The :command:`realpath` command is included in the GNU
               `CoreUtils <https://www.gnu.org/software/coreutils/>`__
               package. Mac OS users can install CoreUtils using
               `MacPorts <https://www.macports.org/>`__ or `Homebrew
               <https://brew.sh/>`__.
