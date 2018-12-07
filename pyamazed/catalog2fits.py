#
# Create a template fits file from a directory of .dat templates
#
usage = """
 usage :
 python catalog2fits.py category source_dir dest_fits

 e.g. :
 python catalog2fits.py galaxy BC03_sdss_tremonti21/galaxy/ galaxy_tmpl.fits
"""

import os
import sys
from astropy.io import fits
from astropy.table import Table

def load_spectrum(name, path):
    wave = []
    flux = []
    with open(path, 'r') as data:
        for line in data:
            if not line or line.startswith('#'):
                continue
            w, f = line.strip().split()
            wave.append(w)
            flux.append(f)
    wave_c = fits.Column(name='WAVE', format='E', array=wave)
    flux_c = fits.Column(name='FLUX', format='E', array=flux)
    hdu = fits.BinTableHDU.from_columns([wave_c, flux_c])
    hdu.name = name
    return hdu

def convert_catalog(category, path, output):
    """
    """
    hdr = fits.Header()
    primary = fits.PrimaryHDU(header=hdr)
    primary.header['CATEGORY'] = category.strip().lower()
    templates = [primary]
    for template in os.scandir(os.path.expanduser(path)):
        if os.path.splitext(template.name)[1] == '.dat':
            print(template.name)
            templates.append(load_spectrum(os.path.splitext(template.name)[0],
                                           template.path))
    hdul = fits.HDUList(templates)
    hdul.writeto(output, overwrite=True)

if __name__ == '__main__':
    convert_catalog(sys.argv[1], sys.argv[2], sys.argv[3])
    try:
        #convert_catalog(sys.argv[1], sys.argv[2], sys.argv[3])
        pass
    except Exception as e:
        print(e)
        print(usage)
