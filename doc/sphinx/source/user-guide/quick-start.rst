***********
Quick Start
***********

Getting started with GYRE is a relatively quick process, involving
just a few steps. For a more in-depth installation guide, which covers
alternative use-cases, refer to the :doc:`installation` chapter.

.. toctree::
   :maxdepth: 2


Install the SDK
===============

First, install the MESA Software Development Kit (SDK) for your
platform (Linux or MacOS) by following the instructions on the `MESA
SDK`_ page.

Build GYRE
==========

Next, download the `GYRE source code`_, and unpack it:

.. substitution-prompt:: bash

   tar xf gyre-|release|.tar.gz

Change into the newly-created ``gyre`` directory, set the ``GYRE_DIR``
environment variable, and build GYRE:

.. substitution-prompt:: bash

   cd gyre-|release|
   export GYRE_DIR=`pwd`
   make -j

(the ``-j`` flags tells ``make`` to use multiple cores, speeding up the build).

Run GYRE
========

Now change into the ``demo`` subdirectory, and launch GYRE to perform
a demonstration calculation:

.. prompt:: bash

   cd demo
   $GYRE_DIR/bin/gyre gyre.in

If all goes well, output will appear on the screen that looks
something like this (edited for brevity):

.. literalinclude:: quick-start/gyre.out
   :language: console
   :end-before: Mode Search
	     
If things go awry, consult the :doc:`troubleshooting`
chapter. Otherwise, proceed to the next chapter where you'll learn
what GYRE actually does.

.. _MESA SDK: mesa-sdk_

.. _GYRE source code: github-tarball_

