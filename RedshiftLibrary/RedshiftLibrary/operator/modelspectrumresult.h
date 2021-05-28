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

  CModelSpectrumResult(const CSpectrum& spc):  m_model(spc) {};
  CModelSpectrumResult(CSpectrum&& spc):  m_model(std::move(spc)) {};

  CModelSpectrumResult() = default;

  CSpectrum& GetSpectrum() { return m_model;};

  void getData(const std::string& name, double **data, int *size) const;
 
private:
  CSpectrum m_model;
};


}

#endif
