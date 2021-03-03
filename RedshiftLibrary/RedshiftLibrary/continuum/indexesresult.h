#ifndef _REDSHIFT_CONTINUUM_INDEXESRESULT_
#define _REDSHIFT_CONTINUUM_INDEXESRESULT_


#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>


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

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;

private:

    Float64 m_StdSpectrum;
    Float64 m_StdContinuum;

};


}

#endif
