// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include <RedshiftLibrary/operator/templatefittingresult.h>
#include <RedshiftLibrary/operator/templatefitting.h>
#include <RedshiftLibrary/extremum/extremum.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <boost/algorithm/string/predicate.hpp>

using namespace NSEpic;

void CTemplateFittingResult::Init(UInt32 n)
{
    ChiSquare.resize( n );
    FitAmplitude.resize( n );
    FitAmplitudeError.resize( n );
    FitAmplitudeSigma.resize( n);
    FitEbmvCoeff.resize( n );
    FitMeiksinIdx.resize( n );
    FitDtM.resize( n );
    FitMtM.resize( n );
    LogPrior.resize( n );
    Redshifts.resize( n );
    Overlap.resize( n );
    Status.resize( n );
    SNR.resize( n );
    ChiSquareIntermediate.resize( n);
    IsmEbmvCoeffIntermediate.resize( n );
    IgmMeiksinIdxIntermediate.resize(n );
}

void CTemplateFittingResult::Init( UInt32 n, Int32 EbmvListSize, Int32 MeiksinListSize)
{
    Init(n);
    InitIntermediate(EbmvListSize, MeiksinListSize);
}

void CTemplateFittingResult::InitIntermediate(Int32 EbmvListSize, Int32 MeiksinListSize)
{ 
    std::vector<TFloat64List> _chi2ListList(EbmvListSize, TFloat64List(MeiksinListSize, DBL_MAX));
    std::vector<TFloat64List> _ismListList(EbmvListSize, TFloat64List(MeiksinListSize, NAN));
    std::vector<TInt32List>   _igmListList(EbmvListSize, TInt32List(MeiksinListSize, -1));

    ChiSquareIntermediate.assign(ChiSquareIntermediate.size(), _chi2ListList);
    IsmEbmvCoeffIntermediate.assign(IsmEbmvCoeffIntermediate.size(), _ismListList);
    IgmMeiksinIdxIntermediate.assign(IgmMeiksinIdxIntermediate.size(), _igmListList);
}

void CTemplateFittingResult::set_at_redshift(const UInt32 i, TFittingIsmIgmResult val)
{
    ChiSquare[i] = val.chiSquare;
    FitAmplitude[i] = val.ampl;
    FitAmplitudeError[i] = val.ampl_err;
    FitAmplitudeSigma[i] = val.ampl_sigma;
    FitEbmvCoeff[i] = val.EbmvCoeff;
    FitMeiksinIdx[i] = val.MeiksinIdx;
    FitDtM[i] = val.sumCross;
    FitMtM[i] = val.sumT;
    LogPrior[i] = val.logprior;
    Overlap[i] = val.overlapRate;
    Status[i] = val.status;

    ChiSquareIntermediate[i] = std::move(val.ChiSquareInterm);
    IsmEbmvCoeffIntermediate[i] = std::move(val.IsmCalzettiCoeffInterm);
    IgmMeiksinIdxIntermediate[i] = std::move(val.IgmMeiksinIdxInterm);
}
