#ifndef _REDSHIFT_LINEMODEL_EXTREMARESULT_
#define _REDSHIFT_LINEMODEL_EXTREMARESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/statistics/pdfcandidatesz.h>

namespace NSEpic
{

class CExtremaResult : public COperatorResult
{

public:

  CExtremaResult() = default;
  CExtremaResult(Int32 n);

  virtual ~CExtremaResult() = default;

  virtual void Resize(Int32 size);
  Int32 size() const;

  TStringList GetIDs() const;
  TFloat64List GetRedshifts() const;
  TFloat64List GetDeltaZs() const;
  TFloat64List GetMerits() const;
  TFloat64List GetValSumProbas() const;

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

  std::string ID(Int32 i) const;
  Float64 Redshift(Int32 i) const;
  Float64 ValProba(Int32 i) const;
  Float64 ValSumProba(Int32 i) const;
  Float64 DeltaZ(Int32 i) const;

  // Extrema results from PDF
  TCandidateZbyRank Candidates;

  //template continuum
  TStringList             FittedTplName;    //Name of the best template fitted for continuum
  TFloat64List            FittedTplAmplitude;     //Amplitude for the best template fitted for continuum
  TFloat64List            FittedTplAmplitudeError;     //Amplitude error for the best template fitted for continuum
  TFloat64List            FittedTplMerit;     //Chisquare for the best template fitted for continuum
  TFloat64List            FittedTplDustCoeff;     //Calzetti dustcoeff for the best template fitted for continuum
  TInt32List        FittedTplMeiksinIdx;    //Meiksin igm index for the best template fitted for continuum
  TFloat64List      FittedTplDtm;    //DTM for the best template fitted for continuum
  TFloat64List      FittedTplMtm;    //MTM for the best template fitted for continuum
  TFloat64List      FittedTplLogPrior;    //log prior for the best template fitted for continuum

protected:
  void SaveJSONbody(std::ostream& stream) const;

};

inline Int32 CExtremaResult::size() const
{
    return Candidates.size();
}

// for compatibility
inline std::string CExtremaResult::ID(Int32 i) const {return Candidates[i].first;}
inline Float64 CExtremaResult::Redshift(Int32 i) const { return Candidates[i].second.Redshift;}
inline Float64 CExtremaResult::ValProba(Int32 i) const { return Candidates[i].second.ValProba;}
inline Float64 CExtremaResult::ValSumProba(Int32 i) const { return Candidates[i].second.ValSumProba;}
inline Float64 CExtremaResult::DeltaZ(Int32 i) const { return Candidates[i].second.Deltaz;}

}

#endif
