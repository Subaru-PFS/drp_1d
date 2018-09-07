from .redshift import *
from astropy.io import fits
import numpy as np
from IPython import embed

class FitsTemplateCatalog(CTemplateCatalog):

    def Load(self, path):
        """Load a template catalog from a fits file"""
        hdul = fits.open(path)

        for spectrum in hdul[1:]:
            name = spectrum.header['EXTNAME']
            print(f'Loading template {name}')

            wavel = spectrum.data.field('WAVE')
            spectralaxis = CSpectrumSpectralAxis(wavel)

            flux = spectrum.data.field('FLUX')
            signal = CSpectrumFluxAxis_withSpectrum(flux)
            template = CTemplate(name, 'galaxy', spectralaxis, signal)
            #embed()
            self.Add(template)
            #template.Save(f'/tmp/dat-{name}')
