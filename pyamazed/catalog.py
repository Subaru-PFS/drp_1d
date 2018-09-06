from .redshift import *
from astropy.io import fits
import numpy as np
from IPython import embed

class FitsTemplateCatalog(CTemplateCatalog):

    def Load(self, path):
        """Load a template catalog from a fits file"""
        hdul = fits.open(path)

        for spectrum in hdul[1:]:
            template = CTemplate()
            wavel = np.array([s[0] for s in spectrum.data])
            spectralaxis = CSpectrumSpectralAxis(wavel)
            flux =  np.array([s[1] for s in spectrum.data])
            #embed()
            signal = CSpectrumFluxAxis_withSpectrum(flux)

            self.Add(template)
