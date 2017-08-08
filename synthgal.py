import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

from .template import Template
from .filters import Filter
import .checkflux as cf


class SynthGal(object):

    def __init__(self, metallicity, haflux, haew, redshift, add_lines, wave_array, temp):

        self.wave_array = wave_array
        # template information
        self.haflux = haflux
        self.haew = haew
        self.redshift = redshift
        self.template = Template(metallicity, haflux, haew, redshift, add_lines, temp)
        #self.mag = self.bandpass_normalization()


    def measure_emission_lines(self):
        """Measures the emission line r"""
#        elements = ['OII', 'OIII', 'Hb', 'OIII', 'OIII', 'NII', 'Ha', 'NII', \
#                    'SII', 'SII', 'SIII', 'SIII']
#        linewaves = [3727., 4364., 4861., 4959., 5007., 6548., 6563., 6583., \
#                     6716., 6730., 9069., 9531.]
        elements = ['OII', 'OIII', 'Hb', 'OIII', 'OIII', 'Ha', \
                    'SII', 'SII', 'SIII', 'SIII']
        linewaves = [3727., 4364., 4861., 4959., 5007., 6563., \
                     6716., 6730., 9069., 9531.]
        strlines = ['%s_%i'%(l,int(w)) for (l,w) in zip(elements,linewaves)]

        # don't want to fit all that craziness
        cond = (self.template.wave > (2500*(1.+self.redshift))) & (self.template.wave < 30000)
        fits,ews,spec = cf.fit_lines(self.template.wave[cond], 
                                     self.template.flux[cond],
                                     linewaves, strlines, self.redshift)

        return zip(strlines, fits, ews)


    def bandpass_normalization(self, filt, fullfilter=True, w1=10400, w2=10600, lam=10500):
        """ """
        # filter used for bandpass normalization
        f = Filter(filt, self.wave_array)
        c = 2.99792e18
        if fullfilter is True:
            print 'using full filter band pass'
            # calculate pivot wavelength and flux density 
            lam,fluxden = f.calc_synphot(self.template.wave, self.template.flux,
                                         self.redshift)
            # calculate magnitude
            mag = -2.5*np.log10(fluxden) - 2.5*np.log10(lam**2) + \
                   2.5*np.log10(c) - 48.6 
        
        else:
            band = (self.template.wave >= w1) & (self.template.wave <= w2)
            fluxden = np.median(self.template.flux[band])
            mag = -2.5*np.log10(fluxden) - 2.5*np.log10(lam**2) + \
                   2.5*np.log10(c) - 48.6 

        return mag
