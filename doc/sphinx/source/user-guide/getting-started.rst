.. _getting_started:

***************
Getting Started
***************

This chapter explains the few simple steps required to get up and
running with GYRE. For a more in-depth installation guide that covers
alternative use-cases, refer to the :ref:`installation` chapter.

.. toctree::
   :maxdepth: 2


Install the SDK
===============

First, install the MESA Software Development Kit (SDK) for your
platform (Linux or MacOS) by following the instructions on the `MESA
SDK <mesa-sdk_>`__ page.

Build GYRE
==========

Next, download the `GYRE source code <github-tarball_>`__, and unpack
it from the command line (these examples assume your using the Bourne
shell):

.. substitution-prompt:: bash

   tar xf gyre-|release|.tar.gz

Change into the newly-created directory, set the ``GYRE_DIR``
environment variable, and build GYRE:

.. substitution-prompt:: bash

   cd gyre-|release|
   export GYRE_DIR=`pwd`
   make -j

(the ``-j`` flags tells ``make`` to use multiple cores, speeding up the build).

Test GYRE
=========

To check that GYRE has compiled correctly and give reasonable results, you can run the calculation test suite via the command

.. substitution-prompt:: bash

   make test

The initial output from the tests should look something like this:

.. literalinclude:: getting-started/make-test.out
   :language: console
   :lines: 1-10

If things go awry, consult the :ref:`troubleshooting`
chapter. Otherwise, proceed to the next chapter where you'll put
together your own GYRE calculation.
