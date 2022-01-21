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
#ifndef _REDSHIFT_SPECTRUM_LOGREBINNING_
#define _REDSHIFT_SPECTRUM_LOGREBINNING_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/processflow/inputcontext.h"

namespace NSEpic
{
class CInputContext;
class CSpectrumLogRebinning 
{

public:
    //applying Rule of zero
    CSpectrumLogRebinning(CInputContext& inputContext);
    std::shared_ptr<CSpectrum>  LoglambdaRebinSpectrum(std::shared_ptr<const CSpectrum> spectrum,
                                                        std::string errorRebinMethod="rebinVariance") const;
    std::shared_ptr<CTemplate>  LoglambdaRebinTemplate(std::shared_ptr<const CTemplate> tpl, 
                                                        TFloat64Range& lambdaRange_tpl, 
                                                        const UInt32  loglambda_count_tpl) const;
    TFloat64Range LogRebinTemplateCatalog(const std::string& category) const;                                                    
    Float64 m_logGridStep; 
    TFloat64Range m_lambdaRange_ref;
private:
    void  SetupRebinning( CSpectrum &spectrum, const TFloat64Range &lambdaRange);
    CSpectrumSpectralAxis  computeTargetLogSpectralAxis(const TFloat64Range & lambdarange, UInt32 gridCount) const;
    UInt32 InferTemplateRebinningSetup(const TFloat64Range& z_range, TFloat64Range& lambdaRange_tpl) const;
    bool CheckTemplateAlignment(const std::shared_ptr<const CTemplate> & tpl, const TFloat64Range& lambdaRange_tpl) const;
    const std::string m_rebinMethod = "lin";

    CInputContext& m_inputContext;
};

}
#endif