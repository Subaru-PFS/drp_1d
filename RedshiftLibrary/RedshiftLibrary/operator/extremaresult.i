
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

  
  
  std::string       FittedTplName ="";    //Name of the best template fitted for continuum
  Float64      FittedTplAmplitude = NAN;     //Amplitude for the best template fitted for continuum
  Float64      FittedTplAmplitudeError= NAN;     //Amplitude error for the best template fitted for continuum
  Float64      FittedTplMerit= NAN;     //Chisquare for the best template fitted for continuum
  Float64      FittedTplEbmvCoeff= NAN;     //Calzetti ebmvcoeff for the best template fitted for continuum
  Int32        FittedTplMeiksinIdx=-1;    //Meiksin igm index for the best template fitted for continuum
  Float64      FittedTplDtm= NAN;    //DTM for the best template fitted for continuum
  Float64      FittedTplMtm= NAN;    //MTM for the best template fitted for continuum
  Float64      FittedTplLogPrior= NAN;    //log prior for the best template fitted for continuum
  Float64      FittedTplSNR= NAN; 
  

};

