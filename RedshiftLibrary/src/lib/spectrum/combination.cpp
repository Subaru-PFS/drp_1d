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
        if(!(spcList[kr]->GetSampleCount()==nSamples))
        {
            sameGrid = false;
            break;
        }

        //check that all spectral axis values are the same
        bool sameSpcAxis = true;
        for(Int32 ks=0; ks<nSamples; ks++)
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

    CSpectrumFluxAxis combinedFlux(nSamples);
    CSpectrumNoiseAxis& combinedNoise = combinedFlux.GetError();
    TAxisSampleList weightSumSq(nSamples);

    for(Int32 kr=0; kr<nRolls; kr++)
    {
        CSpectrumFluxAxis rollFluxAxis = spcList[kr]->GetFluxAxis();
        CSpectrumNoiseAxis& rollNoise = rollFluxAxis.GetError();
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
    
    spcCombined.SetFluxAxis(std::move(combinedFlux));

    return 0;
}
