#include "RedshiftLibrary/operator/modelCombinationcontinuumfittingresult.h"

using namespace NSEpic;

CModelContinuumTplCombinationFittingResult::CModelContinuumTplCombinationFittingResult(
                            Float64 _redshift,                                   
                            Float64 _merit,
                            TFloat64List _amp,
                            TFloat64List _amp_err, 
                            Float64 _ismCoeff,
                            Int32 _igmIndex,
                            Float64 _fitting_snr):
    CModelContinuumFittingResult(_redshift, "undefined", _merit, NAN, NAN, _ismCoeff, _igmIndex, _fitting_snr),
    Amp     ( std::move(_amp)),
    AmpErr  ( std::move(_amp_err))
{}
