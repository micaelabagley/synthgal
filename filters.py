import os
import numpy as np
from astropy.io import fits
from scipy import integrate, interpolate


class Filter(object):
    """A class to define a broadband filter curve.

    The filter transmission curve is read in, normalized, and interpolated
    onto a finer-resolution wavelength array. 
    In :meth:`Filter.calc_synphot`, the spectrum is interpolated 
    to the same wavelength array. The spectrum and the filter used to
    calcuate the broadband photometry therefore have the same wavelength
    resolution.

    Filter transmission curves are stored in the ``configs`` directory
    and are called ``[filter].txt``. Available `HST`/WFC3 filters are:

        - F110W 
        - F160W 
        - F140W

    Args:
        wave_array (float): Wavelength array for interpolation of 
            filters and SED templates
    """

    def __init__(self, filt, wave_array):
        self.wave_array = wave_array
        self.filt = filt
        self.filter_transmission = self.read_filt(self.filt)


    def interp_spectrum(self, wave, flux):
        """Interpolates spectrum to a finer resolution wavelength array.

        Values outside the range of the input wavelength array are set to zero.

        Args:
            wave (float): Input wavelength array
            flux (float): Input flux array

        Returns:
            (float): Flux array interpolated onto higher-res wavelength array
        """
        f = interpolate.interp1d(wave, flux, bounds_error=False, fill_value=0.)
        return f(self.wave_array)


    def read_filt(self, filt):
        """Reads in the filter transmission curve. 

        :func:`read_filt` normalizes the transmission curve and calls 
        :func:`interp_spectrum` to interpolate it to a finer resolution 
        wavelength array.

        Args:
            filt (str): Filter ID

        Returns:
            transmission (float): Interpolated filter transmission curve
        """
        d = np.genfromtxt(os.path.join('configs', filt+'.txt'))
        wave = d[:,0]
        flux = d[:,1]
        # normalize the transmission curve
        norm = integrate.simps(flux, wave)
        flux = flux / norm
        # interpolate transmission curve to a finer wavelength resolution 
        transmission = self.interp_spectrum(wave, flux)
        return transmission


    def calc_synphot(self, wave, flux, redshift):
        """Calculates the flux density of a spectrum in a broad band filter.

        The mean flux density in a broad passband, :math:`P(\lambda)`, 
        is defined as:

        .. math::
        
           f_{\lambda}(P) = \\frac{\int P_{\lambda} \\ f_{\lambda} \\ \lambda \\ \mathrm{ d}\lambda}{\int P_{\lambda} \\ \lambda \\ \mathrm{ d}\lambda}
        
        The pivot wavelength of the filter is then:
        
        .. math::
           \lambda_p(P) = \sqrt{\\frac{\int P(\lambda) \lambda \mathrm{ d}\lambda}{\int P(\lambda) \mathrm{ d}\lambda / \lambda}}

        (See the `Synphot manual <http://www.stsci.edu/institute/software_hardware/stsdas/synphot/SynphotManual.pdf>`_)

        Args:
            wave (float): Rest frame wavelength array of spectrum
            flux (float): Flux array of spectrum
            redshift (float): Source redshift

        Returns:
            (tuple): A tuple containing:
                pivot (float): The pivot wavelength of the filter
                flux_density (float): The mean flux density calculated for the broadband filter
        """
#        # redshift spectrum
#        lobs = wave * (1. + redshift)

        # interpolate spectrum to filter's wavelength resolution
        spectrum = self.interp_spectrum(wave, flux)

        func1 = self.filter_transmission * spectrum * self.wave_array
        func2 = self.filter_transmission * self.wave_array
        func3 = self.filter_transmission / self.wave_array

        flux_density = integrate.simps(func1) / integrate.simps(func2)

        # pivot wavelength of filter
        pivot = np.sqrt(integrate.simps(func2) / integrate.simps(func3))

        return pivot, flux_density



