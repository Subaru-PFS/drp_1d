#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/method/classificationsolve.h>
#include <RedshiftLibrary/method/classificationresult.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/processflow/parameterstore.h>

using namespace NSEpic;

CClassificationSolve::CClassificationSolve(TScopeStack &scope,std::string objectType):
  CSolve("classification",scope,objectType)
{
}

std::shared_ptr<CSolveResult> CClassificationSolve::compute(const CInputContext &inputContext,
                                   COperatorResultStore &resultStore,
                                   TScopeStack &scope)    
{

  
  std::shared_ptr<const CPdfSolveResult> galaxyResult=std::shared_ptr<const CPdfSolveResult>(nullptr); //std::make_shared<const CPdfSolveResult>();
  std::shared_ptr<const CPdfSolveResult> starResult=std::shared_ptr<const CPdfSolveResult>(nullptr);//std::make_shared<const CPdfSolveResult>();
  std::shared_ptr<const CPdfSolveResult> qsoResult=std::shared_ptr<const CPdfSolveResult>(nullptr);

  galaxyResult = std::dynamic_pointer_cast<const CPdfSolveResult>(resultStore.GetGlobalResult("galaxy.result").lock());
  if(inputContext.m_ParameterStore->Get<std::string>("enablestellarsolve") == "yes")
    starResult =  std::dynamic_pointer_cast<const CPdfSolveResult>(resultStore.GetGlobalResult("star.result").lock());
  if(inputContext.m_ParameterStore->Get<std::string>("enableqsosolve") == "yes")
    qsoResult =  std::dynamic_pointer_cast<const CPdfSolveResult>(resultStore.GetGlobalResult("qso.result").lock());

    std::shared_ptr<CClassificationResult> classifResult = std::make_shared<CClassificationResult>();
    
    Float64 qsoLogEvidence = -INFINITY ;
    Float64 stellarLogEvidence = -INFINITY;
    Float64 galaxyLogEvidence = -INFINITY;
    Float64 MaxLogEvidence = -INFINITY;

    if(galaxyResult){
        galaxyLogEvidence = galaxyResult->getEvidence();
        Log.LogInfo( "Found galaxy LogEvidence: %e", galaxyLogEvidence);
        if (galaxyLogEvidence > MaxLogEvidence){
            MaxLogEvidence = galaxyLogEvidence;
            typeLabel = "G";
        }
    }

    if(starResult){
    stellarLogEvidence = starResult->getEvidence();
        Log.LogInfo( "Found stellar LogEvidence: %e", stellarLogEvidence);
        if(stellarLogEvidence > MaxLogEvidence)
        {
            MaxLogEvidence = stellarLogEvidence;
            typeLabel = "S";
        }
    }

    if(qsoResult){
        qsoLogEvidence = qsoResult->getEvidence();
        Log.LogInfo( "Found qso LogEvidence: %e", qsoLogEvidence);
        if(qsoLogEvidence>MaxLogEvidence)
        {
            MaxLogEvidence = qsoLogEvidence;
            typeLabel = "Q";
        }
    }
    Log.LogInfo( "Setting object type: %s", typeLabel.c_str());
    // compute  proba 
    Float64 Pgal = 0.;
    Float64 Pstar = 0.;
    Float64 Pqso = 0.;
    Float64 sum = 0.;

    if (galaxyResult){
        Pgal = exp(galaxyLogEvidence - MaxLogEvidence);
        sum += Pgal;
    }
    if (stellarLogEvidence > -INFINITY){
         Pstar = exp(stellarLogEvidence - MaxLogEvidence);
         sum += Pstar;
    }
    if (qsoLogEvidence > -INFINITY){
        Pqso = exp(qsoLogEvidence - MaxLogEvidence);
        sum += Pqso;
    }
    
    if(sum<=0){
        Log.LogError( "%s: Classification G/S/Q, all probabilities undefined.", __func__);
        throw std::runtime_error("Classification G/S/Q, all probabilities undefined");
    }   

    Pgal /= sum;
    Pstar /= sum;
    Pqso /= sum;
    classifResult->SetTypeLabel(typeLabel);
    classifResult->SetG(galaxyLogEvidence, Pgal);
    classifResult->SetS(stellarLogEvidence, Pstar);
    classifResult->SetQ(qsoLogEvidence, Pqso);

    return classifResult;
}
