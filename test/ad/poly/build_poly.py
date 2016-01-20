#!/usr/bin/env ipython
#
# Build polytrope files

import os
import tempfile as tf

# Build a polytrope file

def build_poly (n_poly, Gamma_1, dxi, toler, filename):

    # Create an input file

    fd, infile = tf.mkstemp()

    f = os.fdopen(fd, 'w')

    f.write('''
&poly
	n_d = 0
	n_poly = {0:24.16e}
        Gamma_1 = {1:24.16e}
/

&num
	dxi = {2:24.16e}
	toler = {3:24.16e}
/

&out
	filename = '{4:s}'
/
'''.format(n_poly, Gamma_1, dxi, toler, filename))

    f.close()

    # Run build_poly

    os.system('./build_poly {0:s}'.format(infile))

    # Delete the input file

    print infile

    #os.remove(infile)

#
            
if __name__ == "__main__":

    Gamma_1 = 1.66666666666666667

    build_poly(0.0, Gamma_1, 0.000244949, 1E-10, '0.0/poly.h5')
    build_poly(1.5, Gamma_1, 0.000365375, 1E-10, '1.5/poly.h5')
    build_poly(3.0, Gamma_1, 0.000689685, 1E-10, '3.0/poly.h5')
    build_poly(4.0, Gamma_1, 0.001497155, 1E-10, '4.0/poly.h5')
