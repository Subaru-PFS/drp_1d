#ifndef _REDSHIFT_OPERATOR_SPECTRARESULT_
#define _REDSHIFT_OPERATOR_SPECTRARESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>


using namespace std;
namespace NSEpic
{

class CSpectraFluxResult : public COperatorResult
{

public:

    CSpectraFluxResult();
    virtual ~CSpectraFluxResult();

    CSpectraFluxResult ( UInt32 optio );

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }
  void getData(const std::string& name, double **data, int *size) const;

    TFloat64List   fluxes;
    TFloat64List   wavel;

    UInt32	         m_optio;

};


}

#endif
