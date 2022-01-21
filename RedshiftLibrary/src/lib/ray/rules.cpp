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
#include "RedshiftLibrary/ray/rules.h"
#include "RedshiftLibrary/ray/linetags.h"
#include "RedshiftLibrary/operator/raydetection.h"

using namespace NSEpic;
using namespace std;
#include <fstream>


CRules::CRules(CSpectrum& spc,  CRayCatalog& detectedCatalog, CRayCatalog& restCatalog, TFloat64Range &lambdaRange, Float64 winsize)
{
    m_spc = spc;
    m_DetectedCatalog = detectedCatalog;
    m_RestCatalog = restCatalog;
    m_lambdaRange = lambdaRange;

    m_winsize = winsize;
}

CRules::~CRules()
{

}

Int32 CRules::check(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet)
{
    Int32 intval = 0;

    bool bval;
    bval = checkRule01(z, matchingSolutionSet);
    if(!bval){
        intval += 1;
    }

    bval = checkRule02(z, matchingSolutionSet);
    if(!bval){
        intval += 10;
    }

    bval = checkRule03(z, matchingSolutionSet);
    if(!bval){
        intval += 100;
    }


    return intval;
}

/**
 * This rule applies when only weak rest lines have been identified in the matching solution.
 * In this case, the absence of strong rest lines must be explained by:
 *  - the strong lines being out of the wavelength range_error
 *  - the string lines being unindentified due to high noise
 *
 */
bool CRules::checkRule01(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet){

    // check if only weak lines are in this solution set
    Int32 nStrong=0;
    for( UInt32 i=0; i<matchingSolutionSet.size(); i++ )
    {
        bool found  = matchingSolutionSet[i].RestRay.GetIsStrong();
        if(found==1){
            nStrong++;
        }
    }
    if(nStrong>0){
        return true;
    }

    //check if the absence of strong lines is justified by the wavelength range
    CRayCatalog::TRayVector strongRestRayList = m_RestCatalog.GetFilteredList(CRay::nType_Emission, CRay::nForce_Strong);
    CRayCatalog::TRayVector strongRestRaysInsideLambdaRangeList;
    Int32 ncatalog = strongRestRayList.size();
    for( UInt32 c=0; c<ncatalog; c++ )
    {
        Float64 lambda = strongRestRayList[c].GetPosition()*(1+z);
        if( lambda >= m_lambdaRange.GetBegin() && lambda <= m_lambdaRange.GetEnd() ){
            strongRestRaysInsideLambdaRangeList.push_back(strongRestRayList[c]);
        }
    }
    if(strongRestRaysInsideLambdaRangeList.size()==0){
        return true;
    }

    //check if the absence of the remaining strong lines is explained by noise
    CLineDetection lineDetection;
    //estimate the weak lines SNR
    Float64 maxSnrWeak = 0.0;
    Float64 maxNoiseWeak = 0.0;

    for( UInt32 i=0; i<matchingSolutionSet.size(); i++ )
    {
        Float64 lambda = matchingSolutionSet[i].RestRay.GetPosition()*(1+z);
        TFloat64Range lambdarange(lambda-m_winsize/2.0, lambda+m_winsize/2.0);
        TInt32Range range = m_spc.GetSpectralAxis().GetIndexesAtWaveLengthRange( lambdarange );

        Float64 flux=0.0;
        Float64 noiseWeak=0.0;
        Float64 snrWeak = lineDetection.ComputeFluxes(m_spc, m_winsize, range, TFloat64List(), &flux, &noiseWeak);

        if( maxSnrWeak < snrWeak){
            maxSnrWeak = snrWeak;
        }
        if( maxNoiseWeak < noiseWeak){
            maxNoiseWeak = noiseWeak;
        }
    }
    //estimate the strong lines SNR
    for( UInt32 c=0; c<strongRestRaysInsideLambdaRangeList.size(); c++ )
    {
        Float64 lambda = strongRestRaysInsideLambdaRangeList[c].GetPosition()*(1+z);
        TFloat64Range lambdarange(lambda-m_winsize/2.0, lambda+m_winsize/2.0);
        TInt32Range range = m_spc.GetSpectralAxis().GetIndexesAtWaveLengthRange( lambdarange );

        Float64 flux=0.0;
        Float64 noise=0.0;
        Float64 snrStrong = lineDetection.ComputeFluxes(m_spc, m_winsize, range, TFloat64List(), &flux, &noise);

        if( snrStrong < maxSnrWeak){
            if( noise < maxNoiseWeak ){
                return false;
            }
        }
    }


    return true;
}

/**
 * This rule applies when one of OIII lines have been detected.
 *
 */
bool CRules::checkRule02(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet){

    // check if the OIII doublet is in this solution set
    Int32 founda=0;
    Int32 foundb=0;
    linetags ltags;
    for( UInt32 i=0; i<matchingSolutionSet.size(); i++ )
    {
        std::string name = matchingSolutionSet[i].RestRay.GetName();
        std::size_t foundstra = name.find(ltags.oIIIa_em);
        if (foundstra!=std::string::npos){
            founda++;
        }
        std::size_t foundstrb = name.find(ltags.oIIIb_em);
        if (foundstrb!=std::string::npos){
            // check if OIIIa would be in the wavelength range
            Float64 lambda = getRestRayLambda(ltags.oIIIa_em)*(1+z);
            if( lambda >= (m_lambdaRange.GetBegin()+m_winsize) && (lambda <= m_lambdaRange.GetEnd()-m_winsize) ){
                foundb++;
            }
        }
    }
    if(foundb==1 && founda==0){
        return false;
    }

    return true;
}

/**
 * This rule applies when Hbeta has been detected.
 *
 */
bool CRules::checkRule03(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet){

    // check if the Hbeta doublet is in this solution set
    Int32 foundHbeta=0;
    Int32 foundHalpha=0;
    linetags ltags;

    //Float64 z = CRayMatchingResult::GetMeanRedshiftSolution(matchingSolutionSet);
    for( UInt32 i=0; i<matchingSolutionSet.size(); i++ )
    {
        std::string name = matchingSolutionSet[i].RestRay.GetName();
        std::size_t foundstr = name.find(ltags.hbeta_em);
        if (foundstr!=std::string::npos){
            // check if Halpha would be in the wavelength range
            Float64 lambda = getRestRayLambda(ltags.halpha_em)*(1+z);
            if( lambda >= (m_lambdaRange.GetBegin()+m_winsize) && (lambda <= m_lambdaRange.GetEnd()-m_winsize) ){
                foundHbeta++;
            }
        }
        foundstr = name.find(ltags.halpha_em);
        if (foundstr!=std::string::npos){
            foundHalpha++;
        }
    }
    if(foundHbeta==1 && foundHalpha!=1){
        return false;
    }

    return true;
}

Float64 CRules::getRestRayLambda(std::string nametag){
    CRayCatalog::TRayVector restRayList = m_RestCatalog.GetFilteredList();
    Int32 ncatalog = restRayList.size();
    for( UInt32 c=0; c<ncatalog; c++ )
    {
        std::string name = restRayList[c].GetName();
        std::size_t foundstr = name.find(nametag.c_str());
        if (foundstr!=std::string::npos){
            return restRayList[c].GetPosition();
        }
    }
    return -1.0;
}
