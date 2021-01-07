#include <RedshiftLibrary/method/classificationsolve.h>
#include <RedshiftLibrary/log/log.h>

using namespace NSEpic;

CClassificationSolve::CClassificationSolve(std::string enableStarFitting, std::string enableQsoFitting)
{
    m_enableStarFitting = enableStarFitting;
    m_enableQsoFitting = enableQsoFitting;
}

CClassificationSolve::~CClassificationSolve()
{

}

void CClassificationSolve::Classify(const CDataStore& store, std::shared_ptr<CSolveResult> galaxyResult, std::shared_ptr<CSolveResult> starResult, std::shared_ptr<CSolveResult> qsoResult)    
{
    classifResult = std::shared_ptr<CClassificationResult>( new CClassificationResult() );
    
    Float64 qsoLogEvidence = -INFINITY ;
    Float64 stellarLogEvidence = -INFINITY;
    Float64 galaxyLogEvidence = -INFINITY;
    Float64 MaxLogEvidence = -INFINITY;

    if(galaxyResult){
        galaxyResult->GetEvidenceFromPdf(store, galaxyLogEvidence);
        Log.LogInfo( "Found galaxy LogEvidence: %e", galaxyLogEvidence);
        if (galaxyLogEvidence > MaxLogEvidence){
            MaxLogEvidence = galaxyLogEvidence;
            typeLabel = "G";
        }
    }

    if(m_enableStarFitting=="yes"){
        Int32 retStellarEv = starResult->GetEvidenceFromPdf(store, stellarLogEvidence);
        if(retStellarEv==0)
        {
            Log.LogInfo( "Found stellar LogEvidence: %e", stellarLogEvidence);
            if(stellarLogEvidence > MaxLogEvidence)
            {
                MaxLogEvidence = stellarLogEvidence;
                typeLabel = "S";
            }
        }
    }
    if(m_enableQsoFitting=="yes"){
        Int32 retQsoEv = qsoResult->GetEvidenceFromPdf(store, qsoLogEvidence);
        if(retQsoEv==0)
        {
            Log.LogInfo( "Found qso LogEvidence: %e", qsoLogEvidence);
            if(qsoLogEvidence>MaxLogEvidence)
            {
                MaxLogEvidence = qsoLogEvidence;
                typeLabel = "Q";
            }
        }
    }

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

    return;
}
