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
#ifndef _REDSHIFT_METHOD_LINEMATCHINGSOLVE_
#define _REDSHIFT_METHOD_LINEMATCHINGSOLVE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/method/linematchingsolveresult.h"
#include "RedshiftLibrary/spectrum/template/template.h"

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 * \class CMethodLineMatchingSolve
 * \brief Solver method based on matching peaks to the lines catalogue.
 */
class CMethodLineMatchingSolve
{

public:

    CMethodLineMatchingSolve();
    ~CMethodLineMatchingSolve();

    std::shared_ptr<CLineMatchingSolveResult> Compute(CDataStore& resultStore,
                                                       const CSpectrum& spc, 
                                                       const TFloat64Range& lambdaRange, 
                                                       const TFloat64Range& redshiftsRange, 
                                                       Float64 redshiftStep, 
                                                       const CRayCatalog& restRayCatalog);

    const std::string GetDescription();

private:

    // Peak Detection
    Float64 m_winsize;
    Float64 m_minsize;
    Float64 m_maxsize;
    Float64 m_detectioncut;
    Float64 m_detectionnoiseoffset;
    Float64 m_cut;
    Float64 m_strongcut;
    Float64 m_enlargeRate;

    // Line Matching
    Bool m_disablegaussianfitqualitycheck;
    Bool m_dynamicLinematching;
    Int64 m_minMatchNum;
    Float64 m_tol;

    // Log
    Bool m_bypassDebug; // If True, debug messages are suppressed even if the --verbose flag is passed.
};

}

#endif
