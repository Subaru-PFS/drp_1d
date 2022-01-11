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
