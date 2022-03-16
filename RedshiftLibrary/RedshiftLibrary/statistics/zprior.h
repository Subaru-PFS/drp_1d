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
#ifndef _REDSHIFT_STATISTICS_ZPRIOR_
#define _REDSHIFT_STATISTICS_ZPRIOR_

#include "RedshiftLibrary/common/datatypes.h"

namespace Statistics_zprior { //boost_test_suite
    //all boost_auto_test_case that use private method
    class NormalizePrior_test;
    class GetStrongLinePresenceLogZPrior_test;
    class GetNLinesSNRAboveCutLogZPrior_test;
}
namespace NSEpic
{

/**
 * \ingroup Redshift
 * CZPrior
 */
class CZPrior
{

public:

    CZPrior(bool normalizePrior, const TFloat64List & redshifts): 
        m_normalizePrior(normalizePrior), m_redshifts(redshifts) {};
    CZPrior(): CZPrior(false, {}) {};

    TFloat64List GetConstantLogZPrior(Int32 nredshifts) const;
    TFloat64List GetStrongLinePresenceLogZPrior(const TBoolList & linePresence, const Float64 penalization_factor) const;
    TFloat64List GetNLinesSNRAboveCutLogZPrior(const TInt32List & nlinesAboveSNR, const Float64 penalization_factor) const;
    TFloat64List GetEuclidNhaLogZPrior(const TFloat64List & redshifts, const Float64 aCoeff) const;
    TFloat64List CombineLogZPrior(const TFloat64List & logprior1, const TFloat64List & logprior2) const;

private:
    friend class Statistics_zprior::NormalizePrior_test;
    friend class Statistics_zprior::GetStrongLinePresenceLogZPrior_test;
    friend class Statistics_zprior::GetNLinesSNRAboveCutLogZPrior_test;

    void NormalizePrior(TFloat64List & logzPrior) const;

    const TFloat64List & m_redshifts;
    bool m_normalizePrior;
};


}

#endif
