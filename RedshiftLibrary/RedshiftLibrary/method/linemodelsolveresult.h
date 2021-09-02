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
#ifndef _REDSHIFT_METHOD_LINEMODELSOLVERESULT_
#define _REDSHIFT_METHOD_LINEMODELSOLVERESULT_

#include "RedshiftLibrary/method/solveresult.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/catalog.h"

#include <memory>
#include <vector>


namespace NSEpic
{


/**
 * \ingroup Redshift
 */
class CLineModelSolveResult : public CPdfSolveResult
{

public:

  CLineModelSolveResult(  const TCandidateZ& BestExtremumResult,
                            const std::string & opt_pdfcombination,
                            Float64 evidence );

    virtual ~CLineModelSolveResult();


/*    Bool GetBestRedshift(Float64& redshift,
                         Float64& merit ,
                         Float64 &sigma,
                         Float64 &snrHa,
                         Float64 &lfHa,
                         Float64 &snrOII,
                         Float64 &lfOII) const;*/
/*    Bool GetBestRedshiftLogArea( Float64& redshift, Float64& merit ) const;*/
/*    Bool GetBestRedshiftFromPdf(Float64& redshift,
                                Float64& merit,
                                Float64& sigma,
                                Float64 &snrHa,
                                Float64 &lfHa,
                                Float64 &snrOII,
                                Float64 &lfOII,
                                std::string &modelTplratio,
                                std::string &modelTplContinuum) const;*/
    // Bool GetRedshiftCandidates( const CDataStore& store,  std::vector<Float64>& redshiftcandidates) const;

/*    void preSave(const CDataStore& store);*/
      //Extrema results
  //  std::shared_ptr<const LineModelExtremaResult> ExtremaResult;

private:

    std::string tplratioName="-1";
    std::string tplcontinuumName="-1";
    Float64 sigma;
    Float64 snrHa=-1.0;
    Float64 lfHa=-1.0;
    Float64 snrOII=-1.0;
    Float64 lfOII=-1.0;

};


}

#endif
