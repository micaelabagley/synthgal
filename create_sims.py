import os
import shutil
from itertools import chain, izip_longest
from glob import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table

from .synthgal import SynthGal
import .runaxesim as ras


def create_dir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)


def prep_directories(depth='deep'):
    """Set up directories and config files for simulations and aXeSim.

    Config files appropriate for G102+F110W and G141+F160W.
    """

    # directories for synthgal
    create_dir('%s/templates'%depth)
    create_dir('%s/simoutputs'%depth)

    # directories for axesim
    create_dir('%s/axesim/CONF'%depth)
    create_dir('%s/axesim/DATA'%depth)
    create_dir('%s/axesim/OUTPUT'%depth)
    create_dir('%s/axesim/OUTSIM'%depth)
    create_dir('%s/axesim/SIMDATA'%depth)

    # config files
    confdir = '%s/axesim/CONF'%depth
    simdatadir = '%s/axesim/SIMDATA'%depth
    shutil.copy('configs/F110W.txt', simdatadir)
    shutil.copy('configs/G102.F110W.V4.32_WISP_6.1.conf', confdir)
    shutil.copy('configs/F160W.txt', simdatadir)
    shutil.copy('configs/G141.F160W.V4.32_WISP_6.1.conf', confdir)
    # sensitivity files
    sens = glob('configs/WFC3.IR.G*sens*.fits')
    for s in sens:
        shutil.copy(s, confdir)



def write_mot(table, output):
    """astropy table writer for sextractor format is not implemented yet"""
    m = open(output, 'w')
#    m.write('# 1  NUMBER\n# 2  X_IMAGE\n# 3  Y_IMAGE\n# 4  A_IMAGE\n# 5  B_IMAGE\n# 6  THETA_IMAGE\n# 7  MAG_F1153W\n# 8  MAG_F1537W\n# 9  SPECTEMP\n# 10 Z\n# 11 MODIMAGE\n')
    m.write('# 1  NUMBER\n# 2  X_IMAGE\n# 3  Y_IMAGE\n# 4  A_IMAGE\n# 5  B_IMAGE\n# 6  THETA_IMAGE\n# 7  MAG_F1050W\n# 8  MAG_F1550W\n# 9  SPECTEMP\n# 10 Z\n# 11 MODIMAGE\n')
    for i,row in enumerate(table):
        t = '\t'.join('%i'%x if isinstance(x,int) else '%.3f'%x for x in row)
        m.write('%s\n'%t)
    m.close()


def simulate_fields(simcat, fullfilt=False, depth='deep'):
    # create wavelength array for interpolation of filters and templates
    wave_array = np.arange(1, 60000.)
    # filter normalization 
    if fullfilt is False:
        tpass_flux_G102 = '1040,1060'
        tpass_flux_G141 = '1540,1560'
    else:
        tpass_flux_G102 = 'F110W.txt'
        tpass_flux_G141 = 'F160W.txt'
        #tpass_flux_G141 = 'F140W.txt'

    pars = np.unique(simcat['par'])
    print 'Simulating %i objects in %i fields' % (simcat['par'].shape[0], 
                                                pars.shape[0])
    for par in pars:
        if par <= 75:
            continue
        # directories for templates and outputs
        create_dir('%s/templates/Par%i'%(depth,par))
        create_dir('%s/simoutputs/Par%i'%(depth,par))
        # simulate one field at a time
        w = np.where(simcat['par'] == par)
        mot = Table(data=None, names=('NUMBER','X_IMAGE','Y_IMAGE','A_IMAGE',
                'B_IMAGE','THETA_IMAGE','MAG_F1153W','MAG_F1537W','SPECTEMP',
                'Z','MODIMAGE'), dtype=[int,float,float,float,float,float,
                float,float,int,float,int])
        t = Table(data=None, names=('NUMBER','flux','EW','z','metallicity',
              'OII_flux','OII_EW','OIII_4363_flux','OIII_4363_EW','Hb_flux',
              'Hb_EW','OIII_4959_flux','OIII_4959_EW','OIII_5007_flux',
              'OIII_5007_EW','Ha_flux','Ha_EW','SII_6716_flux','SII_6716_EW',
              'SII_6730_flux','SII_6730_EW','SIII_9069_flux','SIII_9069_EW',
              'SIII_9531_flux','SIII_9531_EW'),
              dtype=[int,float,float,float,'S5',float,float,float,float,float,
              float,float,float,float,float,float,float,float,float,float,
              float,float,float,float,float])

        for row in simcat[w]:
            templatename = '%s/templates/Par%i/object_%i.dat' % (
                                                        depth,par,row['obj'])
            # add_lines = 1 for all sources
            synthgal = SynthGal(row['metallicity'], row['flux'], row['ew_obs'],
                                row['z'], 1, wave_array, templatename)
            # normalization magnitudes 
            mag1 = synthgal.bandpass_normalization('F110W',fullfilter=fullfilt,
                                                   w1=10400,w2=10600,lam=10500)
            mag2 = synthgal.bandpass_normalization('F160W',fullfilter=fullfilt,
                                                   w1=15400,w2=15600,lam=15500)
            # spectral template index is same as obj ID
            stemp = row['obj']
            # templates have already been redshifted
            redshift = 0.
            # no image templates, using a, b, and theta
            mimage = 0
            mot.add_row([row['obj'], row['x'], row['y'], row['a_image'],
                         row['b_image'], row['t_image'], mag1, mag2, stemp,
                         redshift, mimage])
        
            # check normalization of spectrum
            linefluxes = synthgal.measure_emission_lines()
            strline,lineflux,lineew = zip(*linefluxes)
            lineinfo = np.array([item for item in chain.from_iterable(izip_longest(lineflux,lineew)) if item])
            out = np.append([row['obj'], row['flux'], row['ew_obs'], row['z'], 
                             row['metallicity']], lineinfo)
            t.add_row(out)

        # write out MOT
        mot['MAG_F1153W'].format = '{:.3f}'
        mot['MAG_F1537W'].format = '{:.3f}'
        write_mot(mot, '%s/simoutputs/Par%i/Par%i_mot.dat'%(depth,par,par))

        t.write('%s/simoutputs/Par%i/Par%i_input_table.fits'%(depth,par,par), 
                format='fits', overwrite=True)

        # now do all the axesim things
        previous = glob('%s/axesim/OUTSIM/Par%i_G*_slitless.fits'%(depth,par))

        if ((depth == 'shallow') & (len(previous) == 1)) | ((depth == 'deep') & (len(previous) == 2)):
            a = raw_input('Simulation exists. Rerun axesim? [Y/n] ')
            if a.lower() != 'n':
                print 'Using %s for flux normalization' % tpass_flux_G102
                if depth == 'deep':
                    ras.runaxesim('Par%i'%par, 'G102', 'F110W', gexpt=0, 
                                  dexpt=0, tpass_flux=tpass_flux_G102, 
                                  extraction='YES', depth=depth)
                ras.runaxesim('Par%i'%par, 'G141', 'F160W', gexpt=0, dexpt=0,
                              tpass_flux=tpass_flux_G141, extraction='YES',
                              depth=depth)
            else:
               pass
        else:
            if depth == 'deep':
                print 'Using %s for flux normalization' % tpass_flux_G102
                ras.runaxesim('Par%i'%par, 'G102', 'F110W', gexpt=0, dexpt=0,
                              tpass_flux=tpass_flux_G102, extraction='YES',
                              depth=depth)
            print 'Using %s for flux normalization' % tpass_flux_G141
            ras.runaxesim('Par%i'%par, 'G141', 'F160W', gexpt=0, dexpt=0,
                          tpass_flux=tpass_flux_G141, extraction='YES',
                          depth=depth)


def main():
    depth = 'shallow'
    prep_directories(depth=depth)

    # the list of simulated parameters
    simcat = Table.read('params/simcat_%s.fits'%depth, format='fits')
    simulate_fields(simcat, depth=depth)


if __name__ == '__main__':
    main()
