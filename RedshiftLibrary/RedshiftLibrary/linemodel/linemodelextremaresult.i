
class TLineModelResult : public TExtremaResult
{
public:
  TLineModelResult(const TCandidateZ& candz):
    TExtremaResult(candz)
  {
    this->m_type="TLineModelResult";
  }

    TLineModelResult(const CContinuumModelSolution& cms);
    
    void updateFromContinuumModelSolution(const CContinuumModelSolution& cms,bool all);

    void updateFromLineModelSolution(const CLineModelSolution& cms);

  void updateContinuumFromModel(std::shared_ptr<const CLineModelElementList> lmel);
  void updateTplRatioFromModel(std::shared_ptr<const CLineModelElementList> lmel);

  void updateFromModel(std::shared_ptr<CLineModelElementList> lmel,std::shared_ptr<CLineModelResult> lmresult,bool estimateLeastSquareFast,int indx,const TFloat64Range &lambdaRange,int i_2pass);
 
  Float64            MeritContinuum; //extrema merit for continuum

    Float64            mTransposeM;    // extrema model norm
    Float64            CorrScaleMarg;    // extrema scale marg. correction
    Int32              NDof;   //non zero elements in the lambdarange
    Float64            Redshift_lmfit;// z found with lmfit
    Float64            snrHa;
    Float64            lfHa;
    Float64            snrOII;
    Float64            lfOII;

    Float64            NLinesOverThreshold;  
    Float64            LogArea;   // log area for each extrema
    Float64            LogAreaCorrectedExtrema;   // corrected z for each extrema
    Float64            SigmaZ;    // sigmaz for each extrema

    Float64            StrongELSNR;
    std::vector<std::string>            StrongELSNRAboveCut;
    Float64            bic;    // bayesian information criterion for each extrema
    std::vector<CContinuumIndexes::TContinuumIndexList> ContinuumIndexes; //continuum indexes for each extrema
    CMask      OutsideLinesMask;   //Mask with 0 under the lines and 1 anywhere else
    Float64            OutsideLinesSTDFlux;    //STD measured on the spectrum continuum substracted outside lines
    Float64            OutsideLinesSTDError;   //STD measured on the error spectrum outside lines

    //line width
    Float64      Elv;   //emission line width
    Float64      Alv;   //absorption line width
    std::vector<Float64>      GroupsELv;   //per fitting group line width , EL
    std::vector<Float64>      GroupsALv;   //per fitting group line width , AL

    //template continuum (+ base class)
    Float64      FittedTplRedshift;    //Redshift for the best template fitted for continuum
    std::vector<Float64>      FittedTplpCoeffs;    //poly coeffs for the best template fitted for continuum

    //template ratio
    std::string      FittedTplratioName;   //Name of the best template fitted for tplcorr/tplratio
    Float64      FittedTplratioAmplitude;   //amp of the best template fitted for tplcorr/tplratio
    Float64      FittedTplratioDtm;   //dtm of the best template fitted for tplcorr/tplratio
    Float64      FittedTplratioMtm;   //mtm of the best template fitted for tplcorr/tplratio
    Float64      FittedTplratioIsmCoeff;   //IsmCoeff/EBMV of the best template fitted for tplcorr/tplratio


};

