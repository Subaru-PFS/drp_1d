class TTplCombinationResult: public TExtremaResult 
{

  public:

    TTplCombinationResult(const TCandidateZ& candz):TExtremaResult(candz)
    {
      this->m_type = "TTplCombinationResult";
    }

    TFloat64List              FittedTplAmplitudeList;
    TFloat64List              FittedTplAmplitudeErrorList;
    TFloat64List              FittedTplAmplitudeSigmaList;
    std::vector<TFloat64List> FittedTplCovMatrix;

};

