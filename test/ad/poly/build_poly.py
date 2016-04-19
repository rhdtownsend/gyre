#!/usr/bin/env python
#
# Build polytrope files

import os
import tempfile as tf

# Build a polytrope file

def build_poly (n_poly, Delta_d, xi_d, Gamma_1, dxi, toler, filename):

    n_poly_str = ','.join('{0:24.16e}'.format(n) for n in n_poly)
    Delta_d_str = ','.join('{0:24.16e}'.format(d) for d in Delta_d)
    xi_d_str = ','.join('{0:24.16e}'.format(x) for x in xi_d)

    # Create an input file

    fd, infile = tf.mkstemp()

    f = os.fdopen(fd, 'w')

    f.write('''
&poly
	n_d = {0:d}
	n_poly = {1:s}
        Delta_d = {2:s}
        xi_d = {3:s}
        Gamma_1 = {4:24.16e}
/

&num
	dxi = {5:24.16e}
	toler = {6:24.16e}
/

&out
	filename = '{7:s}'
/
'''.format(len(n_poly)-1, n_poly_str, Delta_d_str, xi_d_str,
           Gamma_1, dxi, toler, filename))

    f.close()

    # Run build_poly

    os.system('./build_poly {0:s}'.format(infile))

    # Delete the input file

#    print(infile)

    os.remove(infile)

#
            
if __name__ == "__main__":

    Gamma_1 = 1.66666666666666667

    build_poly([0.0], [], [], Gamma_1, 0.00244949, 1E-10, '0.0/poly.h5')
    build_poly([0.0,0.0], [0.0], [1.], Gamma_1, 0.00244949, 1E-10, '0.0+0.0/poly.h5')
    build_poly([0.0,0.0], [0.0], [1.], Gamma_1, 0.00244949, 1E-10, '0.0+0.0/poly.h5')
    build_poly([1.5], [], [], Gamma_1, 0.00365375, 1E-10, '1.5/poly.h5')
    build_poly([3.0], [], [], Gamma_1, 0.00689685, 1E-10, '3.0/poly.h5')
    build_poly([3.0,3.0], [0.0], [2.], Gamma_1, 0.00689685, 1E-10, '3.0+3.0/poly.h5')
    build_poly([4.0], [], [], Gamma_1, 0.01497155, 1E-10, '4.0/poly.h5')
