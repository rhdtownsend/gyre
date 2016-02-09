import h5py

def read_poly (file) :

    file = h5py.File(file, 'r')

    xi = file['xi'][...]

    Theta = file['Theta'][...]
    dTheta = file['dTheta'][...]

    n_poly = file.attrs['n_poly']
    Gamma_1 = file.attrs['Gamma_1']

    return {'xi': xi, 'Theta': Theta, 'dTheta': dTheta, 
            'n_poly': n_poly, 'Gamma_1': Gamma_1, 
            'xi_1': xi[-1], 'n': len(xi)}
