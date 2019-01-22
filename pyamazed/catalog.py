from .redshift import CSpectrumSpectralAxis, CSpectrumFluxAxis_withSpectrum, \
    CTemplate, CTemplateCatalog
from astropy.io import fits


class FitsTemplateCatalog(CTemplateCatalog):
    def Load(self, path):
        """Load a template catalog from a fits file"""
        hdul = fits.open(path)
        category = hdul[0].header['CATEGORY']
        for spectrum in hdul[1:]:
            name = spectrum.header['EXTNAME']
            print('Loading {} template {}'.format(category, name))

            wavel = spectrum.data.field('WAVE')
            spectralaxis = CSpectrumSpectralAxis(wavel)

            flux = spectrum.data.field('FLUX')
            signal = CSpectrumFluxAxis_withSpectrum(flux)
            template = CTemplate(name, category, spectralaxis, signal)
            self.Add(template)
