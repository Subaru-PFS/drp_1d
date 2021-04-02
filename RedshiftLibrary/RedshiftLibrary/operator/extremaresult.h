#ifndef _REDSHIFT_OPERATOR_EXTREMARESULT_
#define _REDSHIFT_OPERATOR_EXTREMARESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/statistics/pdfcandidatesz.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>

#include <vector>
#include <memory>

namespace NSEpic
{
class CModelSpectrumResult;
class CModelContinuumFittingResult;

class CExtremaResult : public CPdfCandidateszResult
{

public:

  CExtremaResult() = default;
  CExtremaResult(Int32 n);

  virtual ~CExtremaResult() = default;
  //rule of 5 defaults
  CExtremaResult(const CExtremaResult & ) = default;
  CExtremaResult(CExtremaResult && ) = default;
  CExtremaResult & operator=(const CExtremaResult & ) = default;   
  CExtremaResult & operator=(CExtremaResult && ) = default;   

  virtual void Resize(Int32 size);

  void SaveLine(std::ostream& stream ) const {};
  void Save(std::ostream& stream ) const {};

  virtual void SaveJSON(std::ostream& stream ) const;
  virtual void getCandidateData(const int& rank,const std::string& name, Float64& v) const;
  virtual void getCandidateData(const int& rank,const std::string& name, Int32& v) const;
  virtual void getCandidateData(const int& rank,const std::string& name, std::string& v) const;
  virtual void getCandidateData(const int& rank,const std::string& name, double **data, int *size) const;

  virtual void getData(const std::string& name, Int32& v) const;
  virtual void getData(const std::string& name, Float64& v) const;
  virtual void getData(const std::string& name, std::string& v) const;
  virtual void getData(const std::string& name, double **data, int *size) const;

  //template continuum
  TStringList       FittedTplName;    //Name of the best template fitted for continuum
  TFloat64List      FittedTplAmplitude;     //Amplitude for the best template fitted for continuum
  TFloat64List      FittedTplAmplitudeError;     //Amplitude error for the best template fitted for continuum
  TFloat64List      FittedTplMerit;     //Chisquare for the best template fitted for continuum
  TFloat64List      FittedTplEbmvCoeff;     //Calzetti ebmvcoeff for the best template fitted for continuum
  TInt32List        FittedTplMeiksinIdx;    //Meiksin igm index for the best template fitted for continuum
  TFloat64List      FittedTplDtm;    //DTM for the best template fitted for continuum
  TFloat64List      FittedTplMtm;    //MTM for the best template fitted for continuum
  TFloat64List      FittedTplLogPrior;    //log prior for the best template fitted for continuum
  TFloat64List      FittedTplSNR; 

  std::vector<std::shared_ptr<const CModelSpectrumResult>  > m_savedModelSpectrumResults;
  std::vector<std::shared_ptr<const CModelContinuumFittingResult>  > m_savedModelContinuumFittingResults;

protected:
  void SaveJSONbody(std::ostream& stream) const;

};

}

#endif
