import numpy as np

from cora.util import cubicspline

import simpleforecast

# Value of h from given cosmology
h = 0.67610

# Load in all the files of the background quantities and the power spectrum
#
# Note that all distances are in Mpc, so we need to convert everything
# into Mpc / h (including the PS units)
background = np.genfromtxt('../assumptions/background_zlow.dat')
pk = np.genfromtxt('../assumptions/pk0.dat')


# A convenience function for creating interpolation functions
def interpolater(x, y, log=False):
    if log:
        return cubicspline.LogInterpolater(np.dstack((x, y))[0])
    else:
        return cubicspline.Interpolater(np.dstack((x, y))[0])

# Read in the redshift points
za = background[:, 0]

# Read in each quantity, convert the units and create an interpolation function.
distance = interpolater(za, background[:, 1] * h)
growth = interpolater(za, background[:, 2])
Tb = interpolater(za, background[:, 3] * 1e-3)
bias = interpolater(za, background[:, 4])
number_density = interpolater(za, background[:, 5] / h**3)

# Create a power spectrum function
ps = interpolater(pk[:, 0] / h, pk[:, 1] * h**3, log=True)


class ForecastLowZ(simpleforecast.InterferometerBase):
    """Class for CV forecasting at low-z.

    Overrides various methods to use the standard quantities supplied by Anze.
    """

    def signal_power(self, kpar, kperp, z):

        k = (kpar**2 + kperp**2)**0.5
        mu = kpar / k

        f = 1.0

        return Tb(z)**2 * growth(z)**2 * (bias(z) + f * mu**2)**2 * ps(k)

    def proper_distance(self, z):

        return distance(z)

    def shot_noise(self, z):

        return Tb(z)**2 / number_density(z)


experiment_a = ForecastLowZ()

# Set the array size
experiment_a.num_x = 32
experiment_a.num_y = 32
experiment_a.size_x = 10.0
experiment_a.size_y = 10.0

# Set the frequency range
experiment_a.freq_low = 1420.0 / (1.0 + 6.0)  # Redshift 6
experiment_a.freq_high = 1420.0

# Set the receiver temperature and the observing time
experiment_a.T_recv = 50.0
experiment_a.num_year = 5.0  # not given by Anze


if __name__ == '__main__':

    ## debug anze
    
    sn, nmodes, signal, inv_noise, shot_noise=experiment_a.signal_noise_single(800., 1000., debug=True)
    print nmodes.shape
    print experiment_a.kbin[:10]
    print signal[23:27,23:27]*(0.18694*1e3)**2,'signal'

    #print (np.arange(500)*0.01+0.005)[0:10]
    print sn[23:27,23:27]
    stop()
    
    import h5py

    # Generate and write out the forecast for T_recv=50K
    sn_all = experiment_a.signal_noise_comb()
    sn_split = experiment_a.signal_noise_split()

    with h5py.File('sn_lowz_expA_50K.h5') as fh:

        for ii, (z, sn) in enumerate(sn_split):

            dset = fh.create_dataset('sn_band_%i' % ii, data=sn)
            dset.attrs['z'] = z
            dset.attrs['axes'] = ('k_par', 'k_perp')

        dset = fh.create_dataset('sn_all', data=sn_all)
        dset.attrs['axes'] = ('k_par', 'k_perp')

        fh.create_dataset('k_bin', data=experiment_a.kbin)

    # Write out the forecast for T_recv=10K
    experiment_a.T_recv = 10.0
    sn_all = experiment_a.signal_noise_comb()
    sn_split = experiment_a.signal_noise_split()

    with h5py.File('sn_lowz_expA_10K.h5') as fh:

        for ii, (z, sn) in enumerate(sn_split):

            dset = fh.create_dataset('sn_band_%i' % ii, data=sn)
            dset.attrs['z'] = z

        dset = fh.create_dataset('sn_all', data=sn_all)
        dset.attrs['axes'] = ('k_par', 'k_perp')

        fh.create_dataset('k_bin', data=experiment_a.kbin)
