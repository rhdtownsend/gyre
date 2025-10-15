#!/usr/bin/env python
#
# Build polytrope files

import subprocess

# Build a polytrope file

def build_poly(n_poly, Delta_b, z_b, Gamma_1, dz, toler, filename):

    args = ['./build_poly',
            filename,
            '--n_poly='+','.join('{0:.16e}'.format(n) for n in n_poly),
            f'--Gamma_1={Gamma_1:.16e}',
            f'--dz={dz:.16e}',
            f'--toler={toler:.16e}']

    if len(z_b) > 0:
        args += ['--z_b='+','.join('{0:.16e}'.format(z) for z in z_b)]

    if len(Delta_b) > 0:
        args += ['--Delta_b='+','.join('{0:.16e}'.format(d) for d in Delta_b)]

    subprocess.run(args)

#

if __name__ == "__main__":

    Gamma_1 = 1.66666666666666667

    build_poly([0.0], [], [], Gamma_1, 0.00244949, 1E-10, '0.0/poly.h5')
    build_poly([0.0,0.0], [0.0], [1.], Gamma_1, 0.00244949, 1E-10, '0.0+0.0/poly.h5')
    build_poly([1.5], [], [], Gamma_1, 0.00365375, 1E-10, '1.5/poly.h5')
    build_poly([3.0], [], [], Gamma_1, 0.00689685, 1E-10, '3.0/poly.h5')
    build_poly([3.0,3.0], [0.0], [2.], Gamma_1, 0.00689685, 1E-10, '3.0+3.0/poly.h5')
    build_poly([4.0], [], [], Gamma_1, 0.01497155, 1E-10, '4.0/poly.h5')
