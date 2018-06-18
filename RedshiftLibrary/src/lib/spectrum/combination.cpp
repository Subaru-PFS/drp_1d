#include <RedshiftLibrary/spectrum/combination.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/common/mask.h>

#include <math.h>

using namespace NSEpic;

CSpectrumCombination::CSpectrumCombination()
{
}

CSpectrumCombination::~CSpectrumCombination()
{
}

//return -2 if input rolls are not on the same grid
//return -3 if n input rolls < 1
Int32 CSpectrumCombination::Combine( std::vector<std::shared_ptr<CSpectrum>> spcList, CSpectrum& spcCombined)
{
    Int32 nRolls = spcList.size();
    if(nRolls<1)
    {
        return -3;
    }
    Int32 nSamples = spcList[0]->GetSampleCount();
    //check that the rolls are on the same grid
    bool sameGrid = true;
    for(Int32 kr=1; kr<nRolls; kr++)
    {
        //check that the rolls have the same sample counts
        if(!(spcList[0]->GetSampleCount()==spcList[kr]->GetSampleCount()))
        {
            sameGrid = false;
            break;
        }

        //check that all spectral axis values are the same
        bool sameSpcAxis = true;
        for(Int32 ks=1; ks<spcList[0]->GetSampleCount(); ks++)
        {
            if(spcList[0]->GetSpectralAxis()[ks] != spcList[kr]->GetSpectralAxis()[ks])
            sameSpcAxis=false;
            break;
        }
        if(!sameSpcAxis)
        {
            sameGrid = false;
            break;
        }
    }
    if(!sameGrid)
    {
        return -2;
    }

    //spcCombined = std::shared_ptr<CSpectrum>( new CSpectrum(*spcList[0]) );
    //*spcCombined = *spcList[0];
    Float64* weightSumSq = new Float64 [(int)spcCombined.GetSampleCount()]();
    CSpectrumFluxAxis& combinedFluxAxis = spcCombined.GetFluxAxis();
    Float64* combinedFlux = combinedFluxAxis.GetSamples();
    TFloat64List& combinedNoise = combinedFluxAxis.GetError();
    for(Int32 ks=0; ks<nSamples; ks++)
    {
        weightSumSq[ks]=0.0;
        combinedFlux[ks]=0.0;
        combinedNoise[ks]=0.0;
    }
    for(Int32 kr=0; kr<nRolls; kr++)
    {
        CSpectrumFluxAxis rollFluxAxis = spcList[kr]->GetFluxAxis();
        TFloat64List& rollNoise = rollFluxAxis.GetError();
        for(Int32 ks=0; ks<nSamples; ks++)
        {
            Float64 sigma2Inv = 1.0/(rollNoise[ks]*rollNoise[ks]);
            combinedFlux[ks] += rollFluxAxis[ks]*sigma2Inv;
            weightSumSq[ks] += sigma2Inv;
        }
    }
    for(Int32 ks=0; ks<nSamples; ks++)
    {
        combinedFlux[ks] *= 1.0/weightSumSq[ks];
        combinedNoise[ks] = sqrt(1.0/weightSumSq[ks]);
    }

    delete[] weightSumSq;
    return 0;
}
