import os.path
from astropy.io import fits
from pylibamazed.redshift import (TLSFArguments,
                                  TLSFGaussianVarWidthArgs, TLSFGaussianConstantWidthArgs,
                                  TLSFGaussianConstantResolutionArgs, TLSFGaussianNISPVSSPSF201707Args)



LSFParameters = {
    'GaussianConstantWidth': "width",
    'GaussianConstantResolution': "resolution",
    'GaussianNISPVSSPSF201707': "sourcesize",
    "GaussianVariableWidth": "GaussianVariablewidthFileName",
    "GaussianNISPSIM2016": ""
}


TLSFArgumentsCtor = {'GaussianNISPSIM2016': TLSFArguments,
                     'GaussianConstantWidth': TLSFGaussianConstantWidthArgs,
                     'GaussianConstantResolution': TLSFGaussianConstantResolutionArgs,
                     'GaussianNISPVSSPSF201707': TLSFGaussianNISPVSSPSF201707Args,
                     "GaussianVariableWidth": TLSFGaussianVarWidthArgs}





