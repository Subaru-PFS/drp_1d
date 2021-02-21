#ifndef _REDSHIFT_METHOD_LINEMODELSOLVERESULT_
#define _REDSHIFT_METHOD_LINEMODELSOLVERESULT_

#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/linemodel/linemodelextremaresult.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <memory>
#include <vector>


namespace NSEpic
{

class CProcessFlowContext;

/**
 * \ingroup Redshift
 */
class CLineModelSolveResult : public CSolveResult
{

public:

    enum EType
    {
         nType_raw = 1,
         nType_continuumOnly = 2,
         nType_noContinuum = 3,
         nType_all = 4,
    };

    CLineModelSolveResult();
    virtual ~CLineModelSolveResult();

    void Save(std::ostream& stream ) const;
    void SaveLine(std::ostream& stream ) const;
    Bool GetBestRedshift(Float64& redshift,
                         Float64& merit ,
                         Float64 &sigma,
                         Float64 &snrHa,
                         Float64 &lfHa,
                         Float64 &snrOII,
                         Float64 &lfOII) const;
    Bool GetBestRedshiftLogArea( Float64& redshift, Float64& merit ) const;
    Bool GetBestRedshiftFromPdf(Float64& redshift,
                                Float64& merit,
                                Float64& sigma,
                                Float64 &snrHa,
                                Float64 &lfHa,
                                Float64 &snrOII,
                                Float64 &lfOII,
                                std::string &modelTplratio,
                                std::string &modelTplContinuum) const;
    // Bool GetRedshiftCandidates( const CDataStore& store,  std::vector<Float64>& redshiftcandidates) const;

    Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) ;

    void preSave(const CDataStore& store);
    void getData(const std::string& name, Float64& v) const;

    //Extrema results
    std::shared_ptr<CLineModelExtremaResult> ExtremaResult;

private:


    UInt32 m_bestRedshiftMethod = 2;

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
