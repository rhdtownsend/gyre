.. _namelist-input-files:

********************
Namelist Input Files
********************

This chapter describes the various groups that can appear in the input
files read by the :ref:`GYRE frontends <frontends>`. These files are
in Fortran's `namelist format
<https://cyber.dabamos.de/programming/modernfortran/namelists.html>`__,
a simple text-based format containing one or more namelist
groups. Each group begins with the line :nml:literal:`&name` (where
:nml:literal:`name` is the name of the group); a list of option-value
pairs then follows, and the group ends with a slash :nml:literal:`/`.

.. toctree::
   :maxdepth: 2

   input-files/constants-group
   input-files/grid-group
   input-files/model-group
   input-files/mode-group
   input-files/num-group
   input-files/orbit-group
   input-files/osc-group
   input-files/output-groups
   input-files/rot-group
   input-files/scan-group
   input-files/tide-group
