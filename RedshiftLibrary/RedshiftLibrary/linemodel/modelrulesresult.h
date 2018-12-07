#ifndef _REDSHIFT_LINEMODEL_MODELRULESRESULT_
#define _REDSHIFT_LINEMODEL_MODELRULESRESULT_


#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>


namespace NSEpic
{

  /**
   * \ingroup Redshift
   */
class CModelRulesResult : public COperatorResult
{

public:

    CModelRulesResult( TStringList logStrings );
    CModelRulesResult();
    virtual ~CModelRulesResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

private:

    TStringList LogStrings;
};


}

#endif
