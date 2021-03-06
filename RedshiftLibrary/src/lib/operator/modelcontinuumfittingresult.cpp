#include "RedshiftLibrary/operator/modelcontinuumfittingresult.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include "RedshiftLibrary/spectrum/spectrum.h"

using namespace NSEpic;

/**
 * \brief Attributes values to member variables according to arguments.
 **/
CModelContinuumFittingResult::CModelContinuumFittingResult(Float64 _redshift,
                                                            std::string _name,
                                                            Float64 _merit,
                                                            Float64 _amp,
                                                            Float64 _amp_err, 
                                                            Float64 _ismCoeff,
                                                            Int32 _igmIndex,
                                                            Float64 _fitting_snr):
    Redshift( _redshift),
    Merit   ( _merit),
    Amp     ( _amp),
    AmpErr  ( _amp_err),
    Name    ( _name),
    IsmCoeff( _ismCoeff),
    IgmIndex( _igmIndex),

    Fitting_snr( _fitting_snr)
{

}

/**
 * \brief Empty destructor.
 **/
CModelContinuumFittingResult::~CModelContinuumFittingResult()
{
}

