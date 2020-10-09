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
    CModelSpectrumResult(const CSpectrum& spc, 
                        const std::string tplName,
                        const Float64 dustCoeff,
                        const Int32 meiksinIdx,
                        const Float64 amplitude);
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
  void getData(const std::string& name, Int32& v) const;
  void getData(const std::string& name, std::string& v) const;
  void getData(const std::string& name, Float64& v) const;
 
private:
  CSpectrum m_model;
  std::string m_tplName = "-1";
  Float64 m_amplitude = 0.0;
  Float64 m_dustCoeff = -1.0;
  Int32   m_meiksinIdx = -1.0;
};


}

#endif
