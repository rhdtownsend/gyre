# Python class implementing GYRE's tar_fit_t

import numpy as np
import h5py

import gyre_cheb_fit

class TarFit:

    def __init__ (self, m, k, nu_0, cf):

        self.m = m
        self.k = k
        self.nu_0 = nu_0

        self.cf = cf

    @classmethod
    def load (cls, filename, group='/'):

        with h5py.File(filename, 'r') as f:

            m = f[group].attrs['m']
            k = f[group].attrs['k']
            nu_0 = f[group].attrs['nu_0']

            cf = gyre_cheb_fit.ChebFit.load(filename, 'cf')

        return cls(m, k, nu_0, cf)

    def lam (self, nu):

        x = 2.*np.arctan(nu-self.nu_0)/np.pi

        if self.k >= 0:

            lam = self.cf.eval(x)*self.lam_norm_grav(nu)

        else:

            if self.m > 0:

                if nu <= self.nu_0:
                    lam = self.cf.eval(x)*self.lam_norm_ross(nu)
                else:
                    lam = float('nan')

            elif self.m < 0:

                if nu >= self.nu_0:
                    lam = self.cf.eval(x)*self.lam_norm_ross(nu)
                else:
                    lam = float('nan')

            else:

                lam = float('nan')

        return lam

    def lam_norm_grav (self, nu):

        l = np.abs(self.m) + self.k

        return nu**2 + l*(l+1)

    def lam_norm_ross (self, nu):

        if self.k < -1:
            s = -self.k -1
            return float(self.m)**2/(2*s+1)**2
        elif self.k == -1:
            return nu**2
        else:
            raise Exception('Invalid k')
        
