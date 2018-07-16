#!/usr/bin/env python
#

import numpy as np

# Set up orbital periods

P_orbs = np.linspace(27.5, 32.5, 2500)

with open('gyre_bin.in', 'w') as f:

    # Write header

        f.write('''
&constants
/

&model
	model_type = 'EVOL'
	file = 'iori.mesa'
	file_format = 'MESA'
        uniform_rot = .TRUE.
        Omega_rot = 0.
/

&osc
	nonadiabatic = .TRUE.
/

&num
	diff_scheme = 'COLLOC_GL2'
/

&grid
	n_inner = 5
	alpha_osc = 10
	alpha_exp = 2
/
''')

        # Loop over orbital periods

        for P_orb in P_orbs:

            # Write tide namelist group

            f.write('''
&tide
	freq_orb = {:12.5e}
        freq_orb_units = 'CYC_PER_DAY'
	q = 0.580
	e = 0.745
	omega_static = 0.
	k_max = 50
	l_max = 2
/

'''.format(1./P_orb))
