#ifndef _REDSHIFT_LINEMODEL_LINEMODELEXTREMARESULT_
#define _REDSHIFT_LINEMODEL_LINEMODELEXTREMARESULT_

#include <RedshiftLibrary/operator/extremaresult.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/linemodel/linemodelsolution.h>

//#include <RedshiftLibrary/operator/modelspectrumresult.h>
//#include <RedshiftLibrary/linemodel/modelfittingresult.h>
//#include <RedshiftLibrary/operator/modelcontinuumfittingresult.h>
//#include <RedshiftLibrary/linemodel/modelrulesresult.h>
//#include <RedshiftLibrary/operator/spectraFluxResult.h>


namespace NSEpic
{
class CModelSpectrumResult;
class CModelFittingResult;
class CModelContinuumFittingResult;
class CModelRulesResult;
class CSpectraFluxResult;

class CLineModelExtremaResult : public CExtremaResult
{

public:

  CLineModelExtremaResult() = default;
  CLineModelExtremaResult(Int32 n);
  ~CLineModelExtremaResult() = default;

  void Resize(Int32 size);
  
  void Save(std::ostream& stream ) const {}; 
  void SaveJSON( std::ostream& stream ) const;

  void SaveLine(std::ostream& stream ) const {};
   
  void getCandidateData(const int& rank,const std::string& name, Float64& v) const;
  void getCandidateData(const int& rank,const std::string& name, Int32& v) const;
  void getCandidateData(const int& rank,const std::string& name, std::string& v) const;
  void getCandidateData(const int& rank,const std::string& name, double **data, int *size) const;

  void getData(const std::string& name, Int32& v) const;
  void getData(const std::string& name, Float64& v) const;
  void getData(const std::string& name, std::string& v) const;
  void getData(const std::string& name, double **data, int *size) const;

    //Extrema results
    TFloat64List            MeritContinuum; //extrema merit for continuum

    TFloat64List            mTransposeM;    // extrema model norm
    TFloat64List            CorrScaleMarg;    // extrema scale marg. correction
    TInt32List              NDof;   //non zero elements in the lambdarange
    TFloat64List            Redshift_lmfit;// z found with lmfit
    TFloat64List            snrHa;
    TFloat64List            lfHa;
    TFloat64List            snrOII;
    TFloat64List            lfOII;

    std::vector<TFloat64List> ExtendedRedshifts;    // z range around extrema
    TFloat64List            NLinesOverThreshold;  
    TFloat64List            LogArea;   // log area for each extrema
    TFloat64List            LogAreaCorrectedExtrema;   // corrected z for each extrema
    TFloat64List            SigmaZ;    // sigmaz for each extrema

    TFloat64List            StrongELSNR;
    std::vector<std::vector<std::string>>            StrongELSNRAboveCut;
    TFloat64List            bic;    // bayesian information criterion for each extrema
    std::vector<CContinuumIndexes::TContinuumIndexList> ContinuumIndexes; //continuum indexes for each extrema
    std::vector<CMask>      OutsideLinesMask;   //Mask with 0 under the lines and 1 anywhere else
    TFloat64List            OutsideLinesSTDFlux;    //STD measured on the spectrum continuum substracted outside lines
    TFloat64List            OutsideLinesSTDError;   //STD measured on the error spectrum outside lines

    //line width
    TFloat64List      Elv;   //emission line width
    TFloat64List      Alv;   //absorption line width
    std::vector<TFloat64List>      GroupsELv;   //per fitting group line width , EL
    std::vector<TFloat64List>      GroupsALv;   //per fitting group line width , AL

    //template continuum (+ base class)
    TFloat64List      FittedTplRedshift;    //Redshift for the best template fitted for continuum
    std::vector<TFloat64List>      FittedTplpCoeffs;    //poly coeffs for the best template fitted for continuum

    //template ratio
    std::vector<std::string>      FittedTplratioName;   //Name of the best template fitted for tplcorr/tplratio
    TFloat64List      FittedTplratioAmplitude;   //amp of the best template fitted for tplcorr/tplratio
    TFloat64List      FittedTplratioDtm;   //dtm of the best template fitted for tplcorr/tplratio
    TFloat64List      FittedTplratioMtm;   //mtm of the best template fitted for tplcorr/tplratio
    TFloat64List      FittedTplratioIsmCoeff;   //IsmCoeff/EBMV of the best template fitted for tplcorr/tplratio

    mutable std::map<int,TFloat64List> continuumIndexesColorCopy;
    mutable std::map<int,TFloat64List> continuumIndexesBreakCopy;
    
    std::vector<std::shared_ptr<const CModelFittingResult>  > m_savedModelFittingResults;
    std::vector<std::shared_ptr<const CModelRulesResult>  > m_savedModelRulesResults;
    std::vector<std::shared_ptr<const CSpectraFluxResult>  > m_savedModelContinuumSpectrumResults;

};


}

#endif
