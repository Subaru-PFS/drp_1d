#ifndef _REDSHIFT_CONTINUUM_INDEXESRESULT_
#define _REDSHIFT_CONTINUUM_INDEXESRESULT_


#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>


namespace NSEpic
{

  /**
   * \ingroup Redshift
   */
class CContinuumIndexesResult : public COperatorResult
{

public:

    CContinuumIndexesResult();
    virtual ~CContinuumIndexesResult();

    void SetValues(Float64 stdSpc, Float64 std_continuum);

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

private:

    Float64 m_StdSpectrum;
    Float64 m_StdContinuum;

};


}

#endif
