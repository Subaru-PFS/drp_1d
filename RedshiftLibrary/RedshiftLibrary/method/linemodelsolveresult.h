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
