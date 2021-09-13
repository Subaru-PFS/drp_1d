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
#include "RedshiftLibrary/statistics/zprior.h"

#include <float.h>
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;


CZPrior::CZPrior() {}

CZPrior::~CZPrior() {}


// setting cte priors for all redshift values
TFloat64List CZPrior::GetConstantLogZPrior(UInt32 nredshifts)
{
    TFloat64List logzPrior(nredshifts, 0.0);

    return logzPrior;
}

TFloat64List CZPrior::GetStrongLinePresenceLogZPrior(const TBoolList & linePresence,
                                      const Float64 penalization_factor)
{
    Float64 probaPresent = 1.0;
    Float64 probaAbsent = penalization_factor;
    TFloat64List logzPrior(linePresence.size(), probaAbsent);
    Float64 sum = 0.0;
    for (UInt32 kz = 0; kz < linePresence.size(); kz++)
    {
        if (linePresence[kz])
        {
            logzPrior[kz] = probaPresent;
        } else
        {
            logzPrior[kz] = probaAbsent;
        }
        sum += logzPrior[kz];
    }

    if (sum > 0.0)
    {
        for (UInt32 kz = 0; kz < linePresence.size(); kz++)
        {
            logzPrior[kz] /= sum;
        }
    }

    // switch to log
    for (UInt32 kz = 0; kz < linePresence.size(); kz++)
    {
        logzPrior[kz] = log(logzPrior[kz]);
    }

    return logzPrior;
}

TFloat64List CZPrior::GetNLinesSNRAboveCutLogZPrior(const TInt32List &  nlinesAboveSNR,
                                                  const Float64 penalization_factor)
{
    Int32 nz = nlinesAboveSNR.size();
    Int32 nlinesThres = 2;
    Float64 probaPresent = 1.0;
    Float64 probaAbsent = penalization_factor;
    TFloat64List zPrior(nz, probaAbsent);
    Float64 sum = 0.0;
    for (UInt32 kz = 0; kz < nz; kz++)
    {
        if (nlinesAboveSNR[kz] >= nlinesThres)
        {
            zPrior[kz] = probaPresent;
            //Log.LogDetail("ZPrior: Prior: nlinesAboveSNR[kz] >= nlinesThres for kz=%d", kz);
        } else
        {
            zPrior[kz] = probaAbsent;
        }
        sum += zPrior[kz];
    }

    if (sum > 0.0)
    {
        for (UInt32 kz = 0; kz < nz; kz++)
        {
            zPrior[kz] /= sum;
        }
    }

    // switch to log
    TFloat64List logzPrior(nz, probaAbsent);
    for (UInt32 kz = 0; kz < nz; kz++)
    {
        logzPrior[kz] = log(zPrior[kz]);
    }

    return logzPrior;
}

TFloat64List CZPrior::GetEuclidNhaLogZPrior(const TRedshiftList & redshifts, const Float64 aCoeff)
{
    if(aCoeff <= 0.0)
    {
        Log.LogError( "    CZPrior::GetEuclidNhaLogZPrior: problem found aCoeff<=0: aCoeff=%f", aCoeff);
        throw std::runtime_error("    CZPrior::GetEuclidNhaLogZPrior: problem found aCoeff<=0");
    }
    TFloat64List zPrior(redshifts.size(), 0.0);

    Float64 maxP = -DBL_MAX;
    Float64 minP = DBL_MAX;
    for (UInt32 kz = 0; kz < redshifts.size(); kz++)
    {
        Float64 z = redshifts[kz];
        Float64 z2 = z*z;
        Float64 z3 = z2*z;
        Float64 z4 = z3*z;
        Float64 z5 = z4*z;
        Float64 z6 = z5*z;

        //poly reg pozzetti model 1 at FHa=1e-16
        zPrior[kz] = (- 54.7422088727874*z6
                      + 1203.94994364807*z5
                      - 10409.6716744981*z4
                      + 44240.3837462642*z3
                      - 92914.84430357*z2
                      + 79004.76406*z
                      - 2288.98457865 );

        //shape prior at low z, left of the bell
        Bool enable_low_z_flat = false;
        if(enable_low_z_flat && z<0.7204452872044528){
            zPrior[kz]=20367.877916402278;
        }else if(zPrior[kz]<0){
            zPrior[kz]=DBL_MIN;
        }
        //apply strength
        zPrior[kz] = pow(zPrior[kz], aCoeff);

        if (zPrior[kz] > maxP)
        {
            maxP = zPrior[kz];
        }
        if (zPrior[kz] < minP)
        {
            minP = zPrior[kz];
        }
    }

    Log.LogDebug("Pdfz: zPrior: using HalphaZPrior min=%e", minP);
    Log.LogDebug("Pdfz: zPrior: using HalphaZPrior max=%e", maxP);
    Float64 dynamicCut = 1e12;
    if (maxP > 0.0)
    {
        for (UInt32 kz = 0; kz < redshifts.size(); kz++)
        {
            zPrior[kz] /= maxP;
            if (zPrior[kz] < 1. / dynamicCut)
            {
                zPrior[kz] = 1. / dynamicCut;
            }
        }
    }

    Float64 sum = 0.0;
    for (UInt32 kz = 0; kz < redshifts.size(); kz++)
    {
        sum += zPrior[kz];
    }
    if (sum > 0.0)
    {
        for (UInt32 kz = 0; kz < redshifts.size(); kz++)
        {
            zPrior[kz] /= sum;
        }
    }

    // switch to log
    std::vector<Float64> logzPrior(redshifts.size(), 0.0);
    for (UInt32 kz = 0; kz < redshifts.size(); kz++)
    {
        logzPrior[kz] = log(zPrior[kz]);
    }

    return logzPrior;
}

/**
 * @brief CombineLogZPrior
 * returns a vector with log(prior1*prior2) for each z, with normalization so
 * that sumPrior=1
 * @param logprior1
 * @param logprior2
 * @return
 */
TFloat64List CZPrior::CombineLogZPrior(const TFloat64List & logprior1,
                                     const TFloat64List & logprior2)
{
    bool normalizePrior=true;
    TFloat64List logzPriorCombined;
    if (logprior1.size() != logprior2.size())
    {
        return logzPriorCombined;
    }
    Int32 n = logprior1.size();

    Float64 maxi = DBL_MIN;
    for (UInt32 k = 0; k < n; k++)
    {
        Float64 val = logprior1[k] + logprior2[k];
        if (maxi < val)
        {
            maxi = val;
        }
    }
    std::vector<Float64> zPrior_modif(n, 0.0);
    for (UInt32 k = 0; k < n; k++)
    {
        zPrior_modif[k] = logprior1[k] + logprior2[k] - maxi;
    }

    Float64 sumExpModif = 0.0;
    for (UInt32 k = 0; k < n; k++)
    {
        sumExpModif += exp(zPrior_modif[k]);
    }

    Float64 logSum = log(exp(maxi)) * sumExpModif;
    Float64 lognormterm = 0.0;
    if(normalizePrior)
    {
        lognormterm=logSum;
    }

    logzPriorCombined.resize(n);
    for (Int32 k = 0; k < n; k++)
    {
        logzPriorCombined[k] = logprior1[k] + logprior2[k] - lognormterm;
    }

    return logzPriorCombined;
}

