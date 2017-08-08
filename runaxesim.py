import os
import sys
import shutil
from glob import glob
import re
from pyraf import iraf
from iraf import stsdas, slitless, axe, axesim


def sort_numeric(filename):
    """ """
    base = os.path.basename(filename)
    return int(re.search('\d+',base).group())


def create_spec_list(field, depth):
    """ """
    templatelist = glob(os.path.join(depth,'templates',field,'object_*.dat'))
    templatelist = [os.path.abspath(x) for x in templatelist]
    templatelist.sort(key=sort_numeric)
    f = open('%s/simoutputs/%s/%s_spectra.lis'%(depth,field,field), 'w')
    for template in templatelist:
        f.write('%s\n' % template)
    f.close()


def runaxesim(field, grism, dfilt, gexpt=0, dexpt=0, tpass_flux=None, extraction='NO', depth='deep'):
    """ """
    # create the list of template spectra
    create_spec_list(field, depth)

    mot = '%s_mot.dat'%field
    speclist = '%s_spectra.lis'%field
    shutil.copy(os.path.join('%s/simoutputs/%s'%(depth,field),mot), 
                '%s/axesim/DATA'%depth)
    shutil.copy(os.path.join('%s/simoutputs/%s'%(depth,field),speclist), 
                '%s/axesim'%depth)
    cwd = os.getcwd()
    os.chdir('%s/axesim'%depth)
    os.environ['AXE_IMAGE_PATH'] = './DATA'
    os.environ['AXE_OUTPUT_PATH'] = './OUTPUT'
    os.environ['AXE_CONFIG_PATH'] = './CONF'
    os.environ['AXE_SIMDATA_PATH'] = './SIMDATA'
    os.environ['AXE_OUTSIM_PATH'] = './OUTSIM'

    config = '%s.%s.V4.32_WISP_6.1.conf'%(grism,dfilt)
    output = '%s_%s' % (field, grism)

    if tpass_flux is None:
        tpass_flux = '%s.txt'%dfilt
    
    if gexpt == 0:
        iraf.simdata(incat=mot, config=config, extraction=extraction,
                     tpass_direct='%s.txt'%dfilt, inlist_spec=speclist,
                     tpass_flux=tpass_flux, output_root=output, silent=True)
    else:
        iraf.simdata(incat=mot, config=config, extraction=extraction,
                     tpass_direct='%s.txt'%dfilt, inlist_spec=speclist,
                     tpass_flux=tpass_flux, output_root=output,
                     exptime_disp=gexpt, exptime_dir=dexpt, silent=True)

    """
    page 25 of axesim manual

    "There is NO correction applied if the wavelength given in the MOT and
    the passband specified in tpass_flux are dissimilar, as in this example:
        MAG_F1155
        tpass_flux = [1040,1060]nm"
    So tried removing the F160W flux, so there is only 1 flux column
    
    had to check all spectra to make sure they weren't avg of zero in 
    chosen passband. for F110W, use 12000,13000

    can't use F110W with this iteration of template spectra ... 
        redshifted templates do not cover all of the F110W filter passbannd
    """
    os.chdir(cwd)


def main():
    field = sys.argv[1]
    grism = sys.argv[2]
    dexpt = sys.argv[3]
    gexpt = sys.argv[4]

    if grism == 'G102':
        dfilt = 'F110W'
    elif grism == 'G141':
        dfilt = 'F160W'
    else:
        print 'Unknown grism: %s' %grism
        exit()

    runaxesim(field, grism, dfilt, gexpt, dexpt)


if __name__ == '__main__':
    main()
