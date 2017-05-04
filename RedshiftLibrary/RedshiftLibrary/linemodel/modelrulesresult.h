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

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

private:

    TStringList LogStrings;
};


}

#endif
