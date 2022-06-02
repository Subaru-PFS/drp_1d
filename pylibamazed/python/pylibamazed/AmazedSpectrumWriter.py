from pylibamazed.AbstractSpectrumReader import AbstractSpectrumReader
from astropy.io import fits
import numpy as np
import pandas as pd

class AmazedSpectrumWriter(AbstractSpectrumReader):

    def __init__(self, reader):
        self.reader = reader
        
    def write(self, output_path, sort):
        
        hdr = fits.Header()
        
        d = {'Wave': np.array(self.reader.waves[0]).byteswap().newbyteorder() ,
             'Flux': np.array(self.reader.fluxes[0]).byteswap().newbyteorder(),
             'Err': np.array(self.reader.errors[0]).byteswap().newbyteorder()}
        df = pd.DataFrame(data=d)
    
        is_sorted = all(self.reader.waves[0][i] <= self.reader.waves[0][i+1] for i in range(len(self.reader.waves[0]) - 1))
        if not is_sorted and not sort:
    	    raise Exception(" Spectral axis is not sorted. If you wish to sort, rerun the conversion with --sort option")   
        if not is_sorted and sort:
            df = df.sort_values(by=['Wave'])     
    
        primary = fits.PrimaryHDU(header=hdr)
        wave = fits.Column(name='Wave', format='E', array=df.Wave.values)
        flux = fits.Column(name='Flux', format='E', array=df.Flux.values)
        err = fits.Column(name='Err', format='E', array=df.Err.values)

        
        hdu1_header = fits.Header([fits.Card('SCALE', 1, 'Scale on flux and error'),
                                   fits.Card('WFRAME', self.reader.w_frame,
                                             'air/vacuum')])
        hdu = fits.BinTableHDU.from_columns([wave, flux, err], header=hdu1_header)
        hdu.name = "Spectrum"

        if self.reader.lsf_type == "GaussianConstantWidth":
            lsf_header = fits.Header([fits.Card('LSF_TYPE', 'GaussianConstantWidth', 'Type of LSF'),
                                  fits.Card('S_VALUE', self.reader.lsf_data, 'Scalar value of the LSF, 0 by default')])
        elif self.reader.lsf_type == "no_lsf":
            lsf_header = fits.Header([fits.Card('LSF_TYPE', 'no_lsf', 'Type of LSF'),
                                      fits.Card('S_VALUE', 0, 'Scalar value of the LSF, 0 by default')])
        else:
            raise Exception("Unhandled lsf type :" +self.reader.lsf_type)
        
        lsf_wave = fits.Column(name='Wave', format='E', array=[])
        lsf_width = fits.Column(name='Width', format='E', array=[])
        lsf_hdu = fits.BinTableHDU.from_columns([lsf_wave, lsf_width],
                                            header=lsf_header)
        lsf_hdu.name = "LineSpreadFunction"
        
        photometry_hdu = fits.BinTableHDU.from_columns([], header=fits.Header([]))
        photometry_hdu.name = "Photometry"
        hdul = fits.HDUList([primary, hdu, lsf_hdu, photometry_hdu])
        hdul.writeto(output_path)
