# Python class implementing GYRE's tar_fit_t

import numpy as np
import h5py

import gyre_cheb_fit

class TarFit:

    def __init__ (self, m, k, nu_t, cf_neg, cf_ctr, cf_pos):

        self.m = m
        self.k = k
        self.nu_t = nu_t

        if cf_neg is not None:
            self.cf_neg = cf_neg
        if cf_ctr is not None:
            self.cf_ctr = cf_ctr
        if cf_pos is not None:
            self.cf_pos = cf_pos
        
    @classmethod
    def load (cls, filename, group='/'):

        with h5py.File(filename, 'r') as f:

            m = f[group].attrs['m']
            k = f[group].attrs['k']
            nu_t = f[group].attrs['nu_t']

            if f[group].attrs['df_neg']:
                cf_neg = gyre_cheb_fit.ChebFit.load(filename, 'cf_neg')
            else:
                cf_neg = None

            if f[group].attrs['df_ctr']:
                cf_ctr = gyre_cheb_fit.ChebFit.load(filename, 'cf_ctr')
            else:
                cf_ctr = None

            if f[group].attrs['df_pos']:
                cf_pos = gyre_cheb_fit.ChebFit.load(filename, 'cf_pos')
            else:
                cf_pos = None

        return cls(m, k, nu_t, cf_neg, cf_ctr, cf_pos)

    def lam (self, nu):

        if self.k >= 0:

            if nu <= -self.nu_t:
                lam = self.cf_neg.eval(self.nu_t/nu)*self.lam_norm_grav_outer(nu)
            elif nu >= self.nu_t:
                lam = self.cf_pos.eval(self.nu_t/nu)*self.lam_norm_grav_outer(nu)
            else:
                lam = self.cf_ctr.eval(nu/self.nu_t)*self.lam_norm_grav_inner(nu)

        else:

            if self.m > 0:

                if nu <= -self.nu_t:
                    if self.k == -1:
                        lam = self.cf_neg.eval(self.nu_t/nu)*self.lam_norm_grav_outer(nu)
                    else:
                        lam = self.cf_neg.eval(self.nu_t/nu)*self.lam_norm_ross(nu)
                else:
                    lam = float('nan')

            elif self.m < 0:

                if nu >= self.nu_t:
                    if self.k == -1:
                        lam = self.cf_pos.eval(self.nu_t/nu)*self.lam_norm_grav_outer(nu)
                    else:
                        lam = self.cf_pos.eval(self.nu_t/nu)*self.lam_norm_ross(nu)
                else:
                    lam = float('nan')

            else:

                lam = float('nan')

        return lam

    def lam_norm_grav_inner (self, nu):

        l = np.abs(self.m) + self.k

        return l*(l+1)

    def lam_norm_grav_outer (self, nu):

        if self.m*nu >= 0.:

            if self.k > 0:
                s = self.k - 1
                return nu**2*(2*s + 1)**2
            else:
                return self.m**2

        else:

            s = self.k + 1
            return nu**2*(2*s + 1)**2

    def lam_norm_ross (self, nu):

        s = -self.k -1

        return float(self.m)**2/(2*s+1)**2

