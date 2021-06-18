class TTplCombinationResult: public TExtremaResult 
{

  public:

    TTplCombinationResult(const TCandidateZ& candz):TExtremaResult(candz)
    {
      this->m_type = "TTplCombinationResult";
    }

    TFloat64List              FittedTplAmplitudeList;
    TFloat64List              FittedTplAmplitudeErrorList;
    std::vector<TFloat64List> FittedTplMtmMatrix;

};

