#ifndef _REDSHIFT_LINEMODEL_MODELRULESRESULT_
#define _REDSHIFT_LINEMODEL_MODELRULESRESULT_


#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>


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
