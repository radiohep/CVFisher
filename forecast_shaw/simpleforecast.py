"""Powerspectrum Fisher forecasting.

Uses the simple formalism of arXiv:0910.5007, which is only valid for compact
interferometers, and single-dish/multi-beam experiments.
"""

import numpy as np

from cora.signal import corr21cm
from cora.util import units


class InterferometerBase(object):
    """Base class for forecasting the IM performance on a compact interferometer.

    Attributes
    ----------
    freq_low, freq_high : float
        Band edges in MHz.
    freq_width : float
        Channel width in MHz.
    num_x, num_y : int
        Number of elements in the EW and NS directions.
    size_x, size_y : float
        Physical size of element in metres.
    f_sky : float
        Fraction of sky observed.
    T_recv : float
        Receiver temperature in Kelvin. This is the temperature for a
        single polarisation channel.
    num_year : float
        Number of year of total observation time.
    kmax : float
        Largest k we want to forecast for in h Mpc^{-1}.
    num_k : int
        Number of linear spaced k-bins in each direction.
    """

    freq_high = 1100.0
    freq_low = 300.0

    freq_width = 0.1

    num_x = 32
    num_y = 32

    size_x = 10.0
    size_y = 10.0

    f_sky = 0.5

    T_recv = 50.0

    num_year = 2.0

    kmax = 5.0
    num_k = 500

    def __init__(self):
        self.cr = corr21cm.Corr21cm()

    def signal_power(self, kpar, kperp, z):
        """21cm signal powerspectrum.

        Parameters
        ----------
        kpar, kperp : np.ndarray
            Wavenumbers to calculate at [h / Mpc]
        z : scalar
            Mean redshift of observations.

        Returns
        -------
        pk : np.ndarray
            The power spectrum at the given k [Mpc^3 / h^3]
        """
        return self.cr.powerspectrum(kpar, kperp, z1=z, z2=z)

    def proper_distance(self, z):
        """Distance to redshift z (in Mpc / h)
        """
        return self.cr.cosmology.proper_distance(z)

    def shot_noise(self, z):
        """Shot noise contribution at redshift z (T_b(z)^2 / N(z))
        """
        return 0.0

    def T_sky(self, freq):
        """SKy temperature at a given frequency (in MHz). Uses the fitting formula Anze suggests.
        """
        T_sky = 2000.0 * (freq / 100.0)**(-2.4) + 2.7
        return T_sky

    def beam_size(self, freq):
        """Beam size at specified frequency."""

        return 3e2 / (self.size_x * self.num_x * freq)

    @property
    def beam_num(self):
        """Number of beams formed at any one time."""

        return (2 * self.num_x - 1) * (2 * self.num_y - 1)

    def V_survey(self, freq_low=None, freq_high=None):
        """Volume of survey in (Mpc/h)^3"""

        d = self.proper_distance

        A_sky = 4 * np.pi * self.f_sky

        freq_low = self.freq_low if freq_low is None else freq_low
        freq_high = self.freq_high if freq_high is None else freq_high

        V = A_sky / 3.0 * (d(freq_to_z(freq_low))**3 - d(freq_to_z(freq_high))**3)

        return V

    def window_x(self, x):
        """Tapering of sensitivity in the UV plane (distance specified in m)."""

        return tri(x / (self.size_x * self.num_x))

    def window_y(self, x):
        """Tapering of sensitivity in the UV plane (distance specified in m)."""

        return tri(x / (self.size_y * self.num_y))

    def noise_power(self, freq):
        """Returns a function to calculate the noise PS at the given freq.
        """

        z = freq_to_z(freq)

        beam_size = self.beam_size(freq)
        A_pix = beam_size**2
        A_survey = self.f_sky * 4 * np.pi

        tau = A_pix / A_survey * units.year * self.num_year * self.beam_num

        # Calculate the comoving size of a frequency bin (at the given freq)
        d = self.proper_distance
        dxf = (d(freq_to_z(freq - self.freq_width)) - d(z))

        # Define the window function in k-space for the parallel and perpendicular directions.
        # Use a sinc function for parallel as it is the FT of a top-hat bin. This is probably a bad choice.
        def window_par(kpar):
            y = kpar * dxf / (2 * np.pi)
            return np.sinc(y) * (np.abs(y) < 1.0)

        # Use a triangle for perpendicular window function, this is sensible as it's the UV plane window
        # function for a compact interferometer where each element is uniformly illuminated
        def window_perp(kperp):
            x = (3e2 * kperp * d(z)) / (freq * 2 * np.pi)

            return (self.window_x(x) * self.window_y(x))**0.5  # This is not correct but is probably good for now

        # Calculate the comoving volume of a single pixel (beam)
        V_pix = A_pix * d(z)**2 * dxf

        # Receiver temperature contribution to instrumental Stokes I
        T_recv_I = self.T_recv / 2**0.5

        return inv_noise_ps_21cm(T_recv_I + self.T_sky(freq), tau, V_pix,
                                 self.freq_width, window_par, window_perp)

    def k_grid(self):
        """Determine the grid in k_par and k_perp to use. Override to use a different gridding.
        """
        ka = np.linspace(0.0, self.kmax, self.num_k + 1)
        dk = np.diff(ka)
        km = 0.5 * (ka[1:] + ka[:-1])

        self._kbin = ka
        self._kpar, self._kperp = km[:, np.newaxis], km[np.newaxis, :]
        self._dkpar, self._dkperp = dk[:, np.newaxis], dk[np.newaxis, :]

    _kpar = None

    @property
    def kpar(self):
        """Parallel wavenumber of bin centre [h / Mpc]
        """
        if self._kpar is None:
            self.k_grid()
        return self._kpar

    _kperp = None

    @property
    def kperp(self):
        """Perpendicular wavenumber of bin centre [h / Mpc]
        """
        if self._kperp is None:
            self.k_grid()
        return self._kperp

    _dkpar = None

    @property
    def dkpar(self):
        """Width of bin in parallel direction [h / Mpc]
        """
        if self._dkpar is None:
            self.k_grid()
        return self._dkpar

    _dkperp = None

    @property
    def dkperp(self):
        """Width of bin in perpendicular direction [h / Mpc]
        """
        if self._dkperp is None:
            self.k_grid()
        return self._dkperp

    @property
    def kbin(self):
        """k bins. Elements kbin[i] and kbin[i+1] are the lower and upper bounds of the bins.
        """
        if self._kbin is None:
            self.k_grid()
        return self._kbin

    def signal_noise_single(self, freq_low, freq_high, debug=False):
        """Calculate the S/N for a single frequency band assumed to have
        negligible evolution.

        Parameters
        ----------
        freq_low, freq_high : float
            Lower and upper frequencies of band.
        debug : boolean
            Return more quantities.

        Returns
        -------
        sn : np.ndarray[:, :]
            S/N for the frequency band (i.e PS / sigma(PS)) for each k-bin.
        """
        # Calculate the central redshift
        fc = 0.5 * (freq_low + freq_high)
        zc = freq_to_z(fc)

        # Calculate the power spectrum amplitude across the grid
        signal = self.signal_power(self.kpar, self.kperp, zc)

        # Calculate the number of modes that are averaged over in each power spectrum bin
        V = self.V_survey(freq_low, freq_high)
        nmodes = n_modes_2d(self.dkpar, self.kperp, self.dkperp, V)

        # Calculate the inverse noise power spectrum
        inv_noise = self.noise_power(fc)(self.kpar, self.kperp)

        # Construct the S/N on each PS-bin
        sn = signal * inv_noise / (1.0 + (signal + self.shot_noise(zc)) * inv_noise) * nmodes**0.5

        if debug:
            return sn, nmodes, signal, inv_noise
        else:
            return sn

    def signal_noise_split(self, nband=None):
        """Calculate the S/N for a set of sub-bands.

        Parameters
        ----------
        nband : int, optional
            The number of bands to split into. If not set, calculates the number of bins
            required such that the PS-samples are roughly independent.

        Returns
        -------
        snlist : list of (float, np.ndarray[:, :]) tuples
            List containing central redshift and S/N array pairs.
        """

        # If the number of bins is not explicitly set, calculate the number that gives a large
        # enough frequency length to leave the power spectrum bins roughly uncorrelated
        if nband is None:
            # Calculate the scale corresponding to the parallel bin spacing
            dx = 2 * np.pi / np.median(self.dkpar)
            fc = 0.5 * (self.freq_high + self.freq_low)

            # Numerically estimate dx/df at the centre of the band
            d = self.proper_distance
            dx_df = d(freq_to_z(fc)) - d(freq_to_z(fc + 1.0))

            # Estimate the frequency spacing corresponding to the largest scale
            # and use this to determine the number of bins
            df = dx / dx_df
            nband = int((self.freq_high - self.freq_low) / df) + 1

            print "Largest scale=%f Mpc/h; equivalent freq length=%f MHz; requires %i bins" % (dx, df, nband)

        # Find the edges of the frequency sub-bands
        fw = (self.freq_high - self.freq_low) / nband
        freq_low = self.freq_low + np.arange(nband) * fw
        freq_high = self.freq_low + (np.arange(nband) + 1) * fw

        zc = freq_to_z(0.5 * (freq_low + freq_high))

        # Calculate the S/N for each sub band
        sn = [self.signal_noise_single(fl, fh) for fl, fh in zip(freq_low, freq_high)]

        return zip(zc, sn)

    def signal_noise_comb(self, nbin=None):
        """Power spectrum S/N of the overall frequency band.

        Parameters
        ----------
        nband : int, optional
            Number of bands used internally.

        Returns
        -------
        sn : np.ndarray[:, :]
            S/N on each power spectrum bin.
        """

        snlist = self.signal_noise_split(nbin)

        zc, sn = zip(*snlist)

        return (np.array(sn)**2).sum(axis=0)**0.5


def inv_noise_ps_21cm(T_sys, tau, V_pix, df, window_par, window_perp):
    """Inverse PS noise error for 21cm experiment.

    Parameters are quite low level.

    Parameters
    ----------
    T_sys : scalar
        System temperature
    tau : scalar
        Integration time per pixel in s.
    V_pix : scalar
        Volume of a real space pixel in (Mpc / h)^3.
    df : scalar
        Frequency channel width in MHz.
    window_par, window_perp : functions
        Window functions in parallel and perp directions.

    Returns
    -------
    psfunc : function(kpar, kperp) -> ps
        A function that calculates the inverse noise powerspectrum at the point in k-space.
    """
    noise = T_sys**2 * (V_pix / (tau * df * 1e6))

    def psfunc(kpar, kperp):
        return (window_par(kpar) * window_perp(kperp))**2 / noise

    return psfunc


def n_modes_2d(dkpar, kperp, dkperp, V_survey):
    """Number of independent modes averaged together in the k-space bin.

    Parameters
    ----------
    dkpar : scalar
        Spacing in k-parallel.
    kperp : scalar
        Position in k-perpendicular.
    dkperp : scalar
        Spacing in k-perpendicular.
    V_survey : scalar
        Volume of the survey.

    Returns
    -------
    nmodes : scalar
    """
    return 2 * np.pi * kperp * dkperp * dkpar * V_survey / (2 * (2 * np.pi)**3.0)


def freq_to_z(freq):
    """Calculate a redshift from a frequency.
    """

    return 1420.0 / freq - 1.0


def tri(x):
    """Unit triangle function.
    """
    return np.where(np.abs(x) < 1.0, 1.0 - np.abs(x), np.zeros_like(x))
