#ifndef _REDSHIFT_LINEMODEL_MODELSPECTRUMRESULT_
#define _REDSHIFT_LINEMODEL_MODELSPECTRUMRESULT_


#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

#include <RedshiftLibrary/spectrum/spectrum.h>

namespace NSEpic
{

  /**
   * \ingroup Redshift
   */
class CModelSpectrumResult : public COperatorResult
{

public:

    CModelSpectrumResult(CSpectrum spc);
    CModelSpectrumResult();
    virtual ~CModelSpectrumResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

private:
    CSpectrum model;

};


}

#endif
