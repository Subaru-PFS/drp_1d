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
#include "RedshiftLibrary/spectrum/combination.h"

#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/common/mask.h"

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
