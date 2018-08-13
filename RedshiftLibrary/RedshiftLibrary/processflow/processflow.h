#ifndef _REDSHIFT_PROCESSFLOW_PROCESSFLOW_
#define _REDSHIFT_PROCESSFLOW_PROCESSFLOW_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/processflow/processflow.h>
#include <RedshiftLibrary/operator/operator.h>


#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>

namespace NSEpic
{

class CProcessFlowContext;
class CTemplate;
class CSpectrum;

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
    Int32 getValueFromRefFile( const char* filePath, std::string spcid, Int32 colID, Float64& zref, Int32 reverseInclusion );
};


}

#endif
