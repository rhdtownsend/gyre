.. _working-with-tags:

.. nml:group:: mode
   :no-target:

*****************
Working With Tags
*****************

This chapter uses a worked example to demonstrate tags --- a simple
yet powerful feature allowing a much greater degree of control over
GYRE calculations.

Example Tag Usage
=================

Consider applying :program:`gyre` to calculate the eigenfrequencies of
a red giant branch (RGB) stellar model. Because non-radial p-modes in
the convective envelope couple with high-order g-modes in the
radiative core, the frequency spacing of the non-radial modes is
*much* smaller than that of the radial modes. In such cases, we
ideally want to use a coarse frequency scan for the radial modes and a
fine frequency scan for the non-radial modes.

The following input file, which is designed to work with the :math:`2\,\Msun` RGB
model in :file:`${GYRE_DIR}/models/mesa/rgb/rgb.mesa`, achieves this
goal using tags:

.. literalinclude:: working-with-tags/gyre.in

Observe that each :nml:group:`mode` namelist group has a
:nml:option:`tag` option. When processing a given group,
:program:`gyre` pairs it up with other namelist groups that match
either of the following criteria:

* The :nml:literal:`tag_list` option is not set;
* The :nml:literal:`tag_list` option is set to a comma-separated list,
  and at least one element of this list matches the :nml:option:`tag`
  option.

In the example given above, the :nml:option:`tag_list <osc.tag_list>`
option of the :nml:group:`osc` namelist group is unset; therefore,
this group is paired with all three :nml:group:`mode` namelist groups,
irrespective of their :nml:option:`tag` options. However, the two
:nml:group:`scan` namelist groups each set their :nml:option:`tag_list
<scan.tag_list>` options. In the first group the :nml:value:`radial`
tag appears, and so this group is paired with the first
:nml:group:`mode` namelist group (i.e., the :math:`\ell=0`
mode). Likewise, in the second group the :nml:value:`non-radial` tag
appears, and so this group is paired with the second and third
:nml:group:`mode` namelist groups (i.e., the :math:`\ell=1` and
:math:`\ell=2` modes).

Tag Rules
=========

In addition to the matching criteria given above, there are a couple
of rules that must be obeyed by tags:

* Tag names can't contain commas (however, they can be otherwise arbitrary);
* If a :nml:group:`mode` namelist group doesn't have a
  :nml:option:`tag` option, then only namelists without a
  :nml:literal:`tag_list` option will be paired with it;
* The :nml:group:`constants`, :nml:group:`model`,
  :nml:group:`ad_output`, :nml:group:`nad_output`, and
  :nml:group:`tide_output` namelist groups don't support tags.
