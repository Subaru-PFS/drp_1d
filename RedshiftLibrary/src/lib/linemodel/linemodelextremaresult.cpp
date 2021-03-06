#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"

#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/modelfittingresult.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
using namespace NSEpic;


// TODO this should be a TExtremaResult constructor, using member initialization list
TLineModelResult::TLineModelResult(const CContinuumModelSolution& cms)
{
   FittedTplName = cms.tplName;
    FittedTplAmplitude = cms.tplAmplitude;
    FittedTplAmplitudeError = cms.tplAmplitudeError;
    FittedTplMerit = cms.tplMerit;
    FittedTplEbmvCoeff = cms.tplEbmvCoeff;
    FittedTplMeiksinIdx = cms.tplMeiksinIdx;
    FittedTplRedshift = cms.tplRedshift;
    FittedTplDtm = cms.tplDtm;
    FittedTplMtm = cms.tplMtm;
    FittedTplLogPrior = cms.tplLogPrior;
    FittedTplpCoeffs = cms.pCoeffs;
  }

void TLineModelResult::updateFromContinuumModelSolution(const CContinuumModelSolution& cms,bool all)
  {
    if (all)
      {
    FittedTplName = cms.tplName;
    FittedTplAmplitude = cms.tplAmplitude;
    FittedTplAmplitudeError = cms.tplAmplitudeError;
    FittedTplMerit = cms.tplMerit;
    FittedTplEbmvCoeff = cms.tplEbmvCoeff;
    FittedTplMeiksinIdx = cms.tplMeiksinIdx;
      }
    FittedTplRedshift = cms.tplRedshift;
    FittedTplDtm = cms.tplDtm;
    FittedTplMtm = cms.tplMtm;
    FittedTplLogPrior = cms.tplLogPrior;
    FittedTplpCoeffs = cms.pCoeffs;
  }

void TLineModelResult::updateFromLineModelSolution(const CLineModelSolution& cms)
  {
    Elv= cms.EmissionVelocity;
    Alv= cms.AbsorptionVelocity;
  }

void TLineModelResult::updateContinuumFromModel(std::shared_ptr<const CLineModelElementList> lmel)
  {
    FittedTplName= lmel->getFitContinuum_tplName();
    FittedTplAmplitude= lmel->getFitContinuum_tplAmplitude();
    FittedTplAmplitudeError= lmel->getFitContinuum_tplAmplitudeError();
    FittedTplMerit= lmel->getFitContinuum_tplMerit();
    FittedTplEbmvCoeff= lmel->getFitContinuum_tplIsmEbmvCoeff();
    FittedTplMeiksinIdx= lmel->getFitContinuum_tplIgmMeiksinIdx();
  }

void TLineModelResult::updateTplRatioFromModel(std::shared_ptr<const CLineModelElementList> lmel)
  {
        FittedTplratioName = lmel->getTplshape_bestTplName();
        FittedTplratioIsmCoeff = lmel->getTplshape_bestTplIsmCoeff();
        FittedTplratioAmplitude = lmel->getTplshape_bestAmplitude();
        FittedTplratioDtm = lmel->getTplshape_bestDtm();
        FittedTplratioMtm = lmel->getTplshape_bestMtm();

  }

void TLineModelResult::updateFromModel(std::shared_ptr<CLineModelElementList> lmel,std::shared_ptr<CLineModelResult> lmresult,bool estimateLeastSquareFast,int idx,const TFloat64Range &lambdaRange, int i_2pass)
  {
    
    // TODO : make all these getters const in CLineModelElementList before uncommenting this
    //LineModelSolutions
    Elv = lmel->GetVelocityEmission();
    Alv = lmel->GetVelocityAbsorption();
        
    if (!estimateLeastSquareFast)
      {
        MeritContinuum =
          lmel->getLeastSquareContinuumMerit(lambdaRange);
      } else
      {
        MeritContinuum =
          lmel->getLeastSquareContinuumMeritFast();
      }

    // store model Ha SNR & Flux
    snrHa =
      lmresult->LineModelSolutions[idx].snrHa;
    lfHa =
      lmresult->LineModelSolutions[idx].lfHa;

    // store model OII SNR & Flux
    snrOII =
      lmresult->LineModelSolutions[idx].snrOII;
    lfOII =
      lmresult->LineModelSolutions[idx].lfOII;

    // store the model norm
    mTransposeM =
      lmel->EstimateMTransposeM(lambdaRange);

    // scale marginalization correction
    Float64 corrScaleMarg = lmel->getScaleMargCorrection(); //
    CorrScaleMarg = corrScaleMarg;

    static Float64 cutThres = 3.0;
    Int32 nValidLines = lmresult->GetNLinesOverCutThreshold(i_2pass, cutThres, cutThres);
    NLinesOverThreshold = nValidLines; // m/Float64(1+nValidLines);
    Float64 cumulStrongELSNR = lmel->getCumulSNRStrongEL(); // getStrongerMultipleELAmpCoeff(); //
    StrongELSNR = cumulStrongELSNR;

    std::vector<std::string> strongELSNRAboveCut = lmel->getLinesAboveSNR(3.5);
    StrongELSNRAboveCut = strongELSNRAboveCut;

    Int32 nddl = lmel->GetNElements(); // get the total number of
    // elements in the model
    nddl = lmresult->LineModelSolutions[idx].nDDL; // override nddl by the actual number of elements in
    // the fitted model
    NDof =
      lmel->GetModelNonZeroElementsNDdl();

    Float64 bic = lmresult->ChiSquare[idx] + nddl * log(lmresult->nSpcSamples); // BIC
    // Float64 aic = m + 2*nddl; //AIC
    bic = bic;
    // lmresult->bic = aic + (2*nddl*(nddl+1) )/(nsamples-nddl-1);
    // //AICc, better when nsamples small

    // compute continuum indexes
    // TODO VB is this useful/necessary now ? if there is a computation it should be done before
    //NB AA commented to avoid adding spectrum to getFromModel arguments
    /*
    CContinuumIndexes continuumIndexes;
    ContinuumIndexes =
      continuumIndexes.getIndexes(spectrum, z);
    */
      
    // save the outsideLinesMask
    OutsideLinesMask =
      lmel->getOutsideLinesMask();

    OutsideLinesSTDFlux = lmel->getOutsideLinesSTD(1, lambdaRange);
    OutsideLinesSTDError = lmel->getOutsideLinesSTD(2, lambdaRange);
    Float64 ratioSTD = -1;
    if(OutsideLinesSTDError>0.0)
      {
        ratioSTD = OutsideLinesSTDFlux/OutsideLinesSTDError;
        Float64 ratio_thres = 1.5;
        if(abs(ratioSTD)>ratio_thres || abs(ratioSTD)<1./ratio_thres)
          {
            Log.LogWarning( "  Operator-Linemodel: STD estimations outside lines do not match: ratio=%e, flux-STD=%e, error-std=%e", ratioSTD, OutsideLinesSTDFlux, OutsideLinesSTDError);
          }else{
          Log.LogInfo( "  Operator-Linemodel: STD estimations outside lines found matching: ratio=%e, flux-STD=%e, error-std=%e", ratioSTD, OutsideLinesSTDFlux, OutsideLinesSTDError);
        }
      }else{
      Log.LogWarning( "  Operator-Linemodel: unable to get STD estimations..." );
    }


  }




std::shared_ptr<const COperatorResult> LineModelExtremaResult::getCandidate(const int& rank,const std::string& dataset) const{
      if (dataset == "model_parameters")  return std::make_shared<const TLineModelResult>(this->m_ranked_candidates[rank].second);
      else if (dataset == "fitted_rays")
	{
	  std::shared_ptr<const COperatorResult> cop =  this->m_savedModelFittingResults[rank];
	  return cop;
	}
      else if (dataset == "model")  return this->m_savedModelSpectrumResults[rank];
      else if (dataset == "continuum")  return this->m_savedModelContinuumSpectrumResults[rank];

      else   throw GlobalException(UNKNOWN_ATTRIBUTE,"Unknown dataset");
    }
    
const std::string& LineModelExtremaResult::getCandidateDatasetType(const std::string& dataset) const {
      if (dataset == "model_parameters")      return this->m_ranked_candidates[0].second.getType();
      else if (dataset == "fitted_rays")  return this->m_savedModelFittingResults[0]->getType();
      else if (dataset == "model")  return this->m_savedModelSpectrumResults[0]->getType();
      else if (dataset == "continuum")  return this->m_savedModelContinuumSpectrumResults[0]->getType();
      else   throw GlobalException(UNKNOWN_ATTRIBUTE,"Unknown dataset");
    }

bool LineModelExtremaResult::HasCandidateDataset(const std::string& dataset) const
{
  return (dataset == "model_parameters" || dataset == "model" ||
	  dataset == "continuum" || dataset == "fitted_rays");
}
