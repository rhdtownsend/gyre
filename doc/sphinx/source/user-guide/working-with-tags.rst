.. _working-with-tags:

*****************
Working With Tags
*****************

Tags are a feature that allow a much greater degree of control over
GYRE calculations. They are relatively simple to use, but can be quite
powerful. The simplest way do demonstrate this is via a worked
example.

Example Tag Usage
=================

Consider calculating the eigenfrequencies of a red giant branch (RGB)
stellar model. Because non-radial p-modes in the convective envelope
couple with high-order g-modes in the radiative core, the frequency
spacing of the non-radial modes is *much* smaller than that of the
radial modes. In such cases, we ideally want to use a coarse frequency
scan for the radial modes and a fine frequency scan for the non-radial
modes.

The following input file, which is designed to work with the :math:`2\,\Msun` RGB
model in :file:`${GYRE_DIR}/models/mesa/rgb/rgb.mesa`, achieves this
goal using tags:

.. literalinclude:: working-with-tags/gyre.in

Observe that each :nml_g:`mode` namelist groups has a :nml_n:`tag`
parameter. When processing a given :nml_g:`mode`, GYRE pairs it up
with other namelist groups that match one of the following criteria:

* The namelist group doesn't have a :nml_n:`tag_list` parameter;
* The namelist does have a :nml_n:`tag_list` parameter, *and* the
  parameter value (a comma-separated list) contains the tag value
  defined in the :nml_g:`mode` group.

In the example given above, the :nml_g:`osc` namelist group doesn't
have a :nml_n:`tag_list` parameter; therefore, it is paired with all
three :nml_g:`mode` namelist groups, irrespective of their
:nml_n:`tag` values. However, the two :nml_g:`scan` namelist groups
each have :nml_n:`tag_list` parameters. In the first group the
:nml_v:`radial` tag appears, and so this group is paired with the
first :nml_g:`mode` namelist group (i.e., the :math:`\ell=0`
mode). Likewise, in the second group the :nml_v:`non-radial` tag
appears, and so this group is paired with the second and third
:nml_g:`mode` namelist groups (i.e., the :math:`\ell=1` and
:math:`\ell=2` modes).

Tag Rules
=========

In addition to the matching criteria given above, there are a couple
of rules that must be obeyed by tags:

* Tag names can't contain commas (however, they can be otherwise arbitrary);
* If a :nml_g:`mode` namelist group doesn't have a :nml_n:`tag`
  parameter, then only namelists without a :nml_n:`tag_list` parameter
  will be paired with it;
* The :nml_g:`constants`, :nml_g:`model`, :nml_g:`ad_output` and
  :nml_g:`nad_output` namelist groups don't support tags.
