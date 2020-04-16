#ifndef _REDSHIFT_PROCESSFLOW_PROCESSFLOW_
#define _REDSHIFT_PROCESSFLOW_PROCESSFLOW_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/processflow/parameterstore.h>
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

    template <class T> void SetMedianContinuum(T& continuum, CSpectrum& spectrumWithoutContinuum, CParameterStore& paramStore);
    void SetDFContinuum(const CSpectrum& spectrum, CSpectrum& spectrumWithoutContinuum, CParameterStore& paramStore);
    void SetRawOrZeroMedian(const CSpectrum& spectrum, CSpectrum& spectrumWithoutContinuum, std::string option);
    void EstimateContinuum(const CSpectrum& spectrum, CSpectrum& spectrumWithoutContinuum, CParameterStore& paramStore, CDataStore& dataStore);
    void ProcessRelevance(const CSpectrum& spectrum, const CSpectrum& spectrumWithoutContinuum, CDataStore& dataStore);
    Bool isPdfValid(CProcessFlowContext &ctx) const;
    Int32 getValueFromRefFile( const char* filePath, std::string spcid, Float64& zref, Int32 reverseInclusion );
};


}

#endif
