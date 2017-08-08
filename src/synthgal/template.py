import os
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
from scipy import integrate
import matplotlib.pyplot as plt


class Template(object):
    """A class to read in BC03 templates.

    BC03 templates are stored in the ``synthsrc/highztemplates`` directory 
    and are called ``[IMF]_[metal]_[SFH].fits``, where 

    - IMF is either `salp` or `chab` for Salpeter and Chabrier initial
      mass functions, respectively; 

    - metal is the metallicity ranging from `m22` to `m62` (solar);

    - SFH is one of `tau0p01`, `tau0p50`, and `tau5p00`, for exponentially
      declining star formation histories with characteristic timescales
      of 0.01, 0.5, and 5 Gyr, respectively.
    """

    def __init__(self, metallicity, flux, ew, redshift, add_lines, output):
        self.tempdir = 'templates'
        self.metallicity = metallicity
        self.tempfile = 'bc2003_hr_%s_chab_cSFR_100Myr.spec'%metallicity
        # 10 Myr age for all templates
        self.age = 1.e7
        # emission lines
        self.haflux = flux
        self.haew = ew
        self.redshift = redshift
        self.line_ratios = fits.getdata('configs/line_ratios.fits')
        self.add_lines = add_lines

        self.wave,self.flux = self.read_template()
 
        self.write_template(output)

    
    def find_thing(self, array, value):
        """Returns the array index of the element closest to value."""
        return np.argmin(np.abs(array - value))


    def read_template(self, plotcheck=False):
        """Reads in the BC03 template with the given metallicity.

        Each template FITS table 
        """
        d = np.genfromtxt(os.path.join('templates', self.tempfile))
        wave = d[:,0]
        flux = d[:,1]
        # normalize the spectrum and scale to desired "real" flux units
        lobs,scaled_flux = self.normalize_template(wave, flux)

        # add emission lines to the scaled template
        if self.add_lines == 1:
            template_flux = self.add_emission_lines(lobs, scaled_flux)

        if plotcheck:
            self.plot_spectrum(lobs, scaled_flux, color='r', lw=2)
            self.plot_spectrum(lobs, template_flux, color='k')

        return lobs, scaled_flux
    

    def normalize_template(self, wave, template_flux):
        """ """
        # find 100A-wide region around Halpha
        w1 = self.find_thing(wave, 6500.)
        w2 = self.find_thing(wave, 6600.)
        # template continuum value in this region
        cont = np.median(template_flux[w1:w2])
        # desired continuum given the flux and ew
        desired_cont = self.haflux / self.haew
        # factor required to scale template to the desired "real" units
        scale = desired_cont / cont

        # redshift template
        lobs = wave * (1. + self.redshift)
        # scale template
        scaled_flux = template_flux * scale
        return lobs, scaled_flux


    def add_emission_lines(self, wave, flux):
        """ 
        Line ratios are with respect to Hbeta. 
        Assume Case B for ne = 1e2 cm^-3 and Te = 10000K: Halpha/Hbeta = 2.86

        """
        halpha_wave = 6563. * (1. + self.redshift)
        hbeta_wave = 4861. * (1. + self.redshift)
        hbeta_flux =  self.haflux / 2.86
        # line ratios will depend on metallicity
        if self.metallicity == 'm62':
            key = 'm52_62_72_ratio'
        else:
            key = '%s_ratio' % self.metallicity
        lines = np.append([self.haflux, hbeta_flux], 
                           self.line_ratios[key]*hbeta_flux)
        line_waves = np.array([x*(1.+self.redshift) for x in 
                               self.line_ratios['lambda']])
        waves = np.append([halpha_wave, hbeta_wave], line_waves)

        continuum_line = lambda x, a, b: a*x + b

        a = 1.
        sig = 3.
        for i,line in enumerate(lines):
            line_center = self.find_thing(wave, waves[i])
            # fill in absorption lines before adding emission
            # go out 5 sig on both sides
            c = ((wave > (wave[line_center]-20*sig)) & \
                 (wave < (wave[line_center]-10*sig))) | \
                ((wave > (wave[line_center]+10*sig)) & \
                 (wave < (wave[line_center]+20*sig)))
            lc = (wave > (wave[line_center]-10*sig)) & \
                 (wave < (wave[line_center]+10*sig))
            # at z > 2.3 or so, redder lines land where spectral templates
            # are lower-resolution and too few elements lie within the region
            # around line center. these lines will not be within the grisms
            if (wave[c].shape[0] < 3) & (wave[line_center] > 20000):
                continue
            cont = np.polyfit(wave[c], flux[c], 1)
            continuum = np.polyval(cont, wave[lc])
            flux[lc] = continuum
            gauss = a * np.exp(-(wave-wave[line_center])**2 / (2.*sig**2))
            # normalize line profile so line flux is total flux in line
            norm = integrate.simps(gauss, wave)
            if norm != 0:
                profile = gauss * (line / norm)
                flux += profile

        return flux 


    def write_template(self, output):
        """ """
        t = Table([self.wave, self.flux], names=('wave','flux'))
        ascii.write(t, output, format='commented_header', overwrite=True)


    def plot_spectrum(self, wave, flux, **kwargs):
        """Plots the spectrum as a check."""
        plt.plot(wave, flux, **kwargs)
        plt.xlim(0, 20000)
        
