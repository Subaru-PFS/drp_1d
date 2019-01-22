#include <RedshiftLibrary/ray/rules.h>
#include <RedshiftLibrary/ray/linetags.h>
#include <RedshiftLibrary/operator/raydetection.h>

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

    Bool bval;
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
Bool CRules::checkRule01(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet){

    // check if only weak lines are in this solution set
    Int32 nStrong=0;
    for( UInt32 i=0; i<matchingSolutionSet.size(); i++ )
    {
        Bool found  = matchingSolutionSet[i].RestRay.GetIsStrong();
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
Bool CRules::checkRule02(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet){

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
Bool CRules::checkRule03(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet){

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
