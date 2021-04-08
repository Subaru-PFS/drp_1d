#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/method/reliabilitysolve.h>
#include <RedshiftLibrary/method/reliabilityresult.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/processflow/parameterstore.h>

using namespace NSEpic;

CReliabilitySolve::CReliabilitySolve(TScopeStack &scope,std::string objectType):
  CSolve("reliability",scope,objectType)
{
}

std::shared_ptr<CSolveResult> CReliabilitySolve::compute(std::shared_ptr<const CInputContext> inputContext,
                                                         std::shared_ptr<COperatorResultStore> resultStore,
                                                         TScopeStack &scope)    
{

  
  std::shared_ptr<const CPdfSolveResult> galaxyResult=std::shared_ptr<const CPdfSolveResult>(nullptr); //std::make_shared<const CPdfSolveResult>();

  galaxyResult = std::dynamic_pointer_cast<const CPdfSolveResult>(resultStore->GetGlobalResult("galaxy.result").lock());

  /*
  std::shared_ptr<const CPdfSolveResult> starResult=std::shared_ptr<const CPdfSolveResult>(nullptr);//std::make_shared<const CPdfSolveResult>();
  std::shared_ptr<const CPdfSolveResult> qsoResult=std::shared_ptr<const CPdfSolveResult>(nullptr);

  if(inputContext->GetParameterStore()->Get<std::string>("enablestellarsolve") == "yes")
    starResult =  std::dynamic_pointer_cast<const CPdfSolveResult>(resultStore->GetGlobalResult("star.result").lock());
  if(inputContext->GetParameterStore()->Get<std::string>("enableqsosolve") == "yes")
    qsoResult =  std::dynamic_pointer_cast<const CPdfSolveResult>(resultStore->GetGlobalResult("qso.result").lock());


Float64 qsoLogEvidence = -INFINITY ;
    Float64 stellarLogEvidence = -INFINITY;


    if(starResult){
    stellarLogEvidence = starResult->getEvidence();
        Log.LogInfo( "Found stellar LogEvidence: %e", stellarLogEvidence);
        if(stellarLogEvidence > MaxLogEvidence)
        {
          merit = starResult->getMerit();
        }
    }

    if(qsoResult){
        qsoLogEvidence = qsoResult->getEvidence();
        Log.LogInfo( "Found qso LogEvidence: %e", qsoLogEvidence);
        if(qsoLogEvidence>MaxLogEvidence)
        {
          merit = qsoResult->getMerit();
        }
    }
  */
    std::shared_ptr<CReliabilityResult> reliabResult = std::make_shared<CReliabilityResult>();
  Float64 MaxLogEvidence = -INFINITY;    
    Float64 galaxyLogEvidence = -INFINITY;
    
    Float64 merit = 0;
    if(galaxyResult){
        galaxyLogEvidence = galaxyResult->getEvidence();
        Log.LogInfo( "Found galaxy LogEvidence: %e", galaxyLogEvidence);
        if (galaxyLogEvidence > MaxLogEvidence){
          merit = galaxyResult->getMerit();
        }
    }


    if (std::isnan(merit)) reliabResult->m_ReliabilityLabel="C6";                                           else {
      int reliability = 6 - floor(merit*6);
      if (reliability == 0) reliability = 1;
      std::ostringstream os;                                                                                  os << "C" << reliability;
      reliabResult->m_ReliabilityLabel=os.str();                                                                 }                                     


    return reliabResult;
}
