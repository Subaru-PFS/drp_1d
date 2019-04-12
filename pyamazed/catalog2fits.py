#
# Create a template fits file from a directory of .dat templates
#
import os
import sys
from astropy.io import fits
import argparse


def load_spectrum(category, name, path):
    wave = []
    flux = []
    with open(path, 'r') as data:
        for line in data:
            if not line or line.startswith('#'):
                continue
            w_data = line.strip().split(maxsplit=3)
            wave.append(float(w_data[0]))
            flux.append(float(w_data[1]))
    wave_c = fits.Column(name='WAVE', format='E', array=wave)
    flux_c = fits.Column(name='FLUX', format='E', array=flux)
    hdr = fits.Header()
    hdr['CATEGORY'] = category
    hdu = fits.BinTableHDU.from_columns([wave_c, flux_c])
    hdu.name = name
    hdu.header['CATEGORY'] = category
    return hdu


def convert_catalog(category, path):
    """
    """
    templates = []
    for template in os.scandir(os.path.expanduser(path)):
        if os.path.splitext(template.name)[1] in ['.dat', '.txt']:
            print(template.name)
            templates.append(load_spectrum(category,
                                           os.path.splitext(template.name)[0],
                                           template.path))
    return templates


def create_catalog(output, galaxy=None, star=None, qso=None):
    hdr = fits.Header()
    primary = fits.PrimaryHDU(header=hdr)
    templates = [primary]
    for category, path in (('galaxy', galaxy),
                           ('star', star),
                           ('qso', qso)):
        if path is not None:
            templates.extend(convert_catalog(category, path))
    hdul = fits.HDUList(templates)
    hdul.writeto(output, overwrite=True)


def main():
    parser = argparse.ArgumentParser(
        description='AMAZED template catalog converter.\n\n'
        'Create a FITS template catalog from .dat templates in '
        'calibration directories.',
        epilog='Example : \n\t{} -g BC03_sdss_tremonti21/galaxy/ '
        '-o galaxy_tmpl.fits'.format(os.path.basename(sys.argv[0])),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--output', '-o', metavar='FILE', type=str,
                        help='Output FITS file.', required=True)
    parser.add_argument('--galaxy', '-g', metavar='DIR', type=str,
                        help='Galaxy templates directory')
    parser.add_argument('--star', '-s', metavar='DIR', type=str,
                        help='Star templates directory')
    parser.add_argument('--qso', '-q', metavar='DIR', type=str,
                        help='QSO templates directory')
    args = parser.parse_args()
    if not (args.galaxy or args.star or args.qso):
        print('At least one catalog must be given', file=sys.stderr)
        return
    create_catalog(args.output, args.galaxy, args.star, args.qso)


if __name__ == '__main__':
    main()
