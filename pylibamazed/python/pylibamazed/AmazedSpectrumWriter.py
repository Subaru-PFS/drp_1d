# ============================================================================
#
# This file is part of: AMAZED
#
# Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
#
# https://www.lam.fr/
#
# This software is a computer program whose purpose is to estimate the
# spectrocopic redshift of astronomical sources (galaxy/quasar/star)
# from there 1D spectrum.
#
# This software is governed by the CeCILL-C license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL-C
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-C license and that you accept its terms.
# ============================================================================

import numpy as np
import pandas as pd
from astropy.io import fits


class AmazedSpectrumWriter:

    def __init__(self, reader):
        self.reader = reader

    def write(self, output_path, sort):

        hdr = fits.Header()

        d = {'Wave': np.array(self.reader.waves.get()).byteswap().newbyteorder(),
             'Flux': np.array(self.reader.fluxes.get()).byteswap().newbyteorder(),
             'Err': np.array(self.reader.errors.get()).byteswap().newbyteorder()}
        df = pd.DataFrame(data=d)

        is_sorted = all(self.reader.waves.get()[i] <= self.reader.waves.get()[
                        i + 1] for i in range(len(self.reader.waves.get()) - 1))
        if not is_sorted and not sort:
            raise Exception(
                " Spectral axis is not sorted. If you wish to sort, rerun the conversion with --sort option")
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
            lsf_header = fits.Header([
                fits.Card('LSF_TYPE', 'GaussianConstantWidth', 'Type of LSF'),
                fits.Card('S_VALUE', self.reader.lsf_data, 'Scalar value of the LSF, 0 by default')
            ])
        elif self.reader.lsf_type == "no_lsf":
            lsf_header = fits.Header([
                fits.Card('LSF_TYPE', 'no_lsf', 'Type of LSF'),
                fits.Card('S_VALUE', 0, 'Scalar value of the LSF, 0 by default')
            ])
        else:
            raise Exception("Unhandled lsf type :" + self.reader.lsf_type)

        lsf_wave = fits.Column(name='Wave', format='E', array=[])
        lsf_width = fits.Column(name='Width', format='E', array=[])
        lsf_hdu = fits.BinTableHDU.from_columns([lsf_wave, lsf_width],
                                                header=lsf_header)
        lsf_hdu.name = "LineSpreadFunction"

        photometry_hdu = fits.BinTableHDU.from_columns([], header=fits.Header([]))
        photometry_hdu.name = "Photometry"
        hdul = fits.HDUList([primary, hdu, lsf_hdu, photometry_hdu])
        hdul.writeto(output_path)
