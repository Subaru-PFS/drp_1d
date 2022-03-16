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
#ifndef _REDSHIFT_METHOD_TEMPLATEFITTINGSOLVERESULT_
#define _REDSHIFT_METHOD_TEMPLATEFITTINGSOLVERESULT_

#include "RedshiftLibrary/method/solveresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/catalog.h"

#include <memory>
#include <vector>
#include <unordered_map>
#include <cmath>

namespace NSEpic
{


/**
 * \ingroup Redshift
 */
class CTemplateFittingSolveResult : public CPdfSolveResult
{

public:
    CTemplateFittingSolveResult(const std::string & scope, 
                                const TCandidateZ& ExtremaResult,
                                const std::string & opt_pdfcombination,
                                Float64 evidence );



  //Extrema results

    const std::string m_scope;
    //std::string m_name;

  std::string m_tplName = "undefined";
  Float64 m_amplitude = NAN;
  Float64 m_amplitudeError = NAN;
  Float64 m_EbmvCoeff = NAN;
  Int32   m_meiksinIdx = -1;

  //Not sure it is necessary here
  Float64   m_fittingSNR = NAN;


/*    std::unordered_map<std::string, std::string> m_scope2name = {
        {"templatefittingsolve",      "TemplateFittingSolve"},
        {"templatefittinglogsolve",   "TemplateFittingLogSolve"},
        {"tplcombinationsolve", "TplcombinationSolve"}
    };*/

    //const EType m_type;


};

}

#endif
