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

    void Save(std::ostream& stream ) const;
    void SaveLine(std::ostream& stream ) const;
   

private:

    TStringList LogStrings;
};


}

#endif
