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

import pytest
import numpy as np
from pylibamazed.redshift import CPhotometricData, CPhotometricBand, CPhotBandCatalog

def test_CPhotometricData():
    names = ("bandH", "bandY")
    flux = (1e-14,2e-15)
    fluxerr = (1e-18,3e-18)
    cp = CPhotometricData(names, flux, fluxerr)

    assert cp.GetNameList() == names
    assert cp.GetFlux(names[0])==flux[0]
    assert cp.GetFluxErr(names[1])==fluxerr[1]


def test_CPhotometricBand():
    lbda = np.arange(12000.,12500.,13.)
    trans = np.full(lbda.size, .7)
    pb = CPhotometricBand(trans,lbda)

    assert np.all(pb.GetTransmission() == trans)
    assert np.all(pb.GetWavelength() == lbda)


def test_CPhotBandCatalog():
    lbda1 = np.arange(12000.,12500.,13.)
    trans1 = np.full(lbda1.size, .7)
    pb1 = CPhotometricBand(trans1,lbda1)

    lbda2 = np.arange(14000.,14300.,13.)
    trans2 = np.full(lbda2.size, .6)
    pb2 = CPhotometricBand(trans2,lbda2)

    names = ("band1", "band2") 
    photcat = CPhotBandCatalog()
    photcat.Add(names[0], pb1)
    photcat.Add(names[1], pb2)

    assert photcat.GetNameList() == names
    pb = photcat[names[0]]

    assert np.all(np.array(pb.GetTransmission())==np.array(pb1.GetTransmission()))
    assert np.all(np.array(pb.GetWavelength())==np.array(pb1.GetWavelength()))
