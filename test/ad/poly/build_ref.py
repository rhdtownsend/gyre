#!/usr/bin/env python
#
# Build reference files for polytropes

import numpy as np
import h5py as h5

# Write a file in GYRE's summary format

def write_file (file, l, n, omega) :

    # Write the file

    h5.get_config().complex_names = ('re', 'im')

    f = h5.File(file, 'w')

    f.attrs['label'] = np.string_(' '*256)

    f.create_dataset('l', data=l)
    f.create_dataset('n_pg', data=n)
    f.create_dataset('omega', data=omega)

    f.close()

    
# Build a reference file using analytical frequencies (from [Pek1938])

def build_ref_0 (out_file, gamma) :

    # Create the data

    l_min = 0
    l_max = 3
    n_min = 1
    n_max = 10

    N = (l_max-l_min+1)*(n_max-n_min+1)

    l = np.empty(N, dtype=np.int32)
    n = np.empty(N, dtype=np.int32)
    omega = np.empty(N, dtype=np.complex128)

    j = 0
    
    for l_ in range(l_min, l_max+1) :
        for n_ in range(n_min, n_max+1) :

            l[j] = l_
            n[j] = n_

            k = 2*(n_-1)

            D = 0.5*(-4. + 0.5*gamma*(k*(k+5+2*l_) + 6 + 4*l_))

            beta = D + np.sqrt(D**2 + l_*(l_+1))

            omega[j] = np.sqrt(beta)

            j += 1

    # Write the output file

    write_file(out_file, l, n, omega)

    
# Build a reference file using tabulated frequency data

def build_ref (in_file, out_file, nu_ref, l_0) :

    # Read the input file

    d = np.genfromtxt(in_file, dtype=[('i4'), ('f8'), ('f8'), ('f8'), ('f8')])

    # Process the data

    l_min = l_0
    l_max = l_0 + 3

    N = (l_max-l_min+1)*len(d)

    l = np.empty(N, dtype=np.int32)
    n = np.empty(N, dtype=np.int32)
    omega = np.empty(N, dtype=np.complex128)

    j = 0
    
    for l_ in range(l_min, l_max+1) :
        for i in range(len(d)):
            l[j] = l_
            n[j] = d['f0'][i]
            omega[j] = d['f{:d}'.format(l_-l_min+1)][i]/nu_ref
            j += 1

    # Write the output file

    write_file(out_file, l, n, omega)


#

if __name__ == '__main__':

    gamma = 5./3.

    build_ref_0('0.0/ref/summary.h5', gamma)

    nu_ref = 99.855377

    build_ref('data/poly-freqs-1.5.dat', '1.5/ref/summary.h5', nu_ref, 0)
    build_ref('data/poly-freqs-3.0.dat', '3.0/ref/summary.h5', nu_ref, 0)
    build_ref('data/poly-freqs-3.0-g.dat', '3.0-g/ref/summary.h5', nu_ref, 1)
    build_ref('data/poly-freqs-4.0.dat', '4.0/ref/summary.h5', nu_ref, 0)
