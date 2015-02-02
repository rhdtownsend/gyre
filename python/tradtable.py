import h5py
import numpy as np
import numpy.polynomial.chebyshev as cheby

import trad

class TradTable :

    """Table of traditional-approximation data"""

    def __init__ (self, filename) :

        """Initialize from an HDF5 file"""

        f = h5py.File(filename, 'r')

        self.l_max = f.attrs['l_max']
        self.trad = np.empty([self.l_max+1,2*self.l_max+1], dtype=object)

        for l in range(0, self.l_max+1) :
            for m in range(-l, l+1) :

                tf_group_name = 'tf({0:d},{1:d})'.format(l, m)

                def read_cheby (cb_group_name) :
                    c = f[cb_group_name]['c'][...]
                    x_a = f[cb_group_name].attrs['x_a']
                    x_b = f[cb_group_name].attrs['x_b']
                    c[0] *= 0.5
                    return cheby.Chebyshev(c, domain=[x_a,x_b])

                cb_pos = read_cheby('{0:s}/cb_pos'.format(tf_group_name))
                cb_neg = read_cheby('{0:s}/cb_neg'.format(tf_group_name))
                cb_ctr = read_cheby('{0:s}/cb_ctr'.format(tf_group_name))

                assert f[tf_group_name].attrs['m'] == m
                k = f[tf_group_name].attrs['k']

                self.trad[l,l+m] = trad.Trad(m, k, cb_pos, cb_neg, cb_ctr)
                
        f.close()

    def lamb (self, l, m, nu) :

        """Laplace tidal equation eigenvalue"""

        return self.trad[l,l+m].lamb(nu)

    def l_e (self, l, m, nu) :

        """Effective harmonic degree"""

        return self.trad[l,l+m].l_e(nu)







        
