import h5py

def read_poly (file) :

    file = h5py.File(file, 'r')

    z = file['z'][...]

    theta = file['theta'][...]
    dtheta = file['dtheta'][...]

    n_poly = file.attrs['n_poly']
    Gamma_1 = file.attrs['Gamma_1']

    return {'z': z, 'theta': theta, 'dtheta': dtheta, 
            'n_poly': n_poly, 'Gamma_1': Gamma_1, 
            'z_s': xi[-1], 'n': len(xi)}
