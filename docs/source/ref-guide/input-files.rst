.. _namelist-input-files:

********************
Namelist Input Files
********************

This chapter describes the various groups that can appear in the input
files read by the :ref:`GYRE frontends <frontends>`. These files are
in Fortran's `namelist format
<https://cyber.dabamos.de/programming/modernfortran/namelists.html>`__,
a simple text-based format containing one or more namelist
groups. Each group begins with the line :nml_g:`name` (where ``name``
is the name of the group); a list of parameter-value pairs then
follows, and the group ends with a slash ``/``.

.. toctree::
   :maxdepth: 2

   input-files/constants
   input-files/grid-params
   input-files/model-params
   input-files/mode-params
   input-files/num-params
   input-files/orbit-params
   input-files/osc-params
   input-files/output-params
   input-files/rot-params
   input-files/scan-params
   input-files/tidal-params
