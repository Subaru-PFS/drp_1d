#ifndef _REDSHIFT_LINEMODEL_MODELSPECTRUMRESULT_
#define _REDSHIFT_LINEMODEL_MODELSPECTRUMRESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/common/ref.h>

namespace NSEpic
{

class CModelSpectrumResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CModelSpectrumResult )

public:

    CModelSpectrumResult(CSpectrum spc);
    CModelSpectrumResult();
    virtual ~CModelSpectrumResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;

private:
    CSpectrum model;

};


}

#endif
