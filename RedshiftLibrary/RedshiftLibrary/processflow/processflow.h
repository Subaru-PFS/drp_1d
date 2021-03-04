#ifndef _REDSHIFT_PROCESSFLOW_PROCESSFLOW_
#define _REDSHIFT_PROCESSFLOW_PROCESSFLOW_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/operator.h>

namespace NSEpic
{

class CProcessFlowContext;

/**
 * \ingroup Redshift
 */
class CProcessFlow
{

public:

    CProcessFlow();
    ~CProcessFlow();

    void Process( CProcessFlowContext& ctx );

private:

    Bool isPdfValid(CProcessFlowContext &ctx) const;
    Int32 getValueFromRefFile( const char* filePath, std::string spcid, Float64& zref, Int32 reverseInclusion );
};


}

#endif
