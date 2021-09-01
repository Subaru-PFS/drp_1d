// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/method/solveresult.h"
#include "RedshiftLibrary/method/classificationsolve.h"
#include "RedshiftLibrary/method/classificationresult.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/processflow/parameterstore.h"

using namespace NSEpic;

CClassificationSolve::CClassificationSolve(TScopeStack &scope,std::string objectType):
  CSolve("classification",scope,objectType)
{
}

std::shared_ptr<CSolveResult> CClassificationSolve::compute(std::shared_ptr<const CInputContext> inputContext,
                                                            std::shared_ptr<COperatorResultStore> resultStore,
                                                            TScopeStack &scope)    
{

  
  std::shared_ptr<const CPdfSolveResult> galaxyResult=std::shared_ptr<const CPdfSolveResult>(nullptr); //std::make_shared<const CPdfSolveResult>();
  std::shared_ptr<const CPdfSolveResult> starResult=std::shared_ptr<const CPdfSolveResult>(nullptr);//std::make_shared<const CPdfSolveResult>();
  std::shared_ptr<const CPdfSolveResult> qsoResult=std::shared_ptr<const CPdfSolveResult>(nullptr);

  galaxyResult = std::dynamic_pointer_cast<const CPdfSolveResult>(resultStore->GetSolveResult("galaxy"));
  if(inputContext->GetParameterStore()->Get<std::string>("enablestellarsolve") == "yes")
    starResult =  std::dynamic_pointer_cast<const CPdfSolveResult>(resultStore->GetSolveResult("star"));
  if(inputContext->GetParameterStore()->Get<std::string>("enableqsosolve") == "yes")
    qsoResult =  std::dynamic_pointer_cast<const CPdfSolveResult>(resultStore->GetSolveResult("qso"));

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
