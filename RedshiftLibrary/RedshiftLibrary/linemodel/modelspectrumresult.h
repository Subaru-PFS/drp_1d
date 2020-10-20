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

    CModelSpectrumResult(const CSpectrum& spc);
    CModelSpectrumResult();
    virtual ~CModelSpectrumResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }
    CSpectrum& GetSpectrum();

  void getData(const std::string& name, double **data, int *size) const;
  

  
private:
    CSpectrum m_model;

};


}

#endif
