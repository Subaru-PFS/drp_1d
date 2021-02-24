
class TExtremaResult : public TCandidateZ {

public:
  TExtremaResult()
    {
    };
  TExtremaResult(const TExtremaResult&) = default;
  TExtremaResult(TExtremaResult&&) = default;
  TExtremaResult& operator=(const TExtremaResult&) = default;
  TExtremaResult& operator=(TExtremaResult&&) = default;
  virtual ~TExtremaResult() = default;
 TExtremaResult(const TCandidateZ& candz):
  TCandidateZ(candz)
    {
      this->m_type = "TExtremaResult";
    }
  
  std::string       FittedTplName;    //Name of the best template fitted for continuum
  Float64      FittedTplAmplitude;     //Amplitude for the best template fitted for continuum
  Float64      FittedTplAmplitudeError;     //Amplitude error for the best template fitted for continuum
  Float64      FittedTplMerit;     //Chisquare for the best template fitted for continuum
  Float64      FittedTplEbmvCoeff;     //Calzetti ebmvcoeff for the best template fitted for continuum
  Int32        FittedTplMeiksinIdx;    //Meiksin igm index for the best template fitted for continuum
  Float64      FittedTplDtm;    //DTM for the best template fitted for continuum
  Float64      FittedTplMtm;    //MTM for the best template fitted for continuum
  Float64      FittedTplLogPrior;    //log prior for the best template fitted for continuum
  Float64      FittedTplSNR; 
  

};

