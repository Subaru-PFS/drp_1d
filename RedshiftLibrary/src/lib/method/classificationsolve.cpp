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
    
    Float64 qsoEvidence = 0 ;
    Float64 stellarEvidence = 0;
    Float64 galaxyEvidence = 0;
    
    if(galaxyResult){
        galaxyResult->GetEvidenceFromPdf(store, galaxyEvidence);
        Log.LogInfo( "Found galaxy evidence: %e", galaxyEvidence);
        typeLabel = "G";
    }

    if(m_enableStarFitting=="yes"){
        Int32 retStellarEv = starResult->GetEvidenceFromPdf(store, stellarEvidence);
        if(retStellarEv==0)
        {
            Log.LogInfo( "Found stellar evidence: %e", stellarEvidence);
            if(stellarEvidence>galaxyEvidence)
            {
                typeLabel = "S";
            }
        }
    }
    if(m_enableQsoFitting=="yes"){
        Int32 retQsoEv = qsoResult->GetEvidenceFromPdf(store, qsoEvidence);
        if(retQsoEv==0)
        {
            Log.LogInfo( "Found qso evidence: %e", qsoEvidence);
            if(qsoEvidence>galaxyEvidence && qsoEvidence>stellarEvidence)
            {
                typeLabel = "Q";
            }
        }
    }
    classifResult->SetTypeLabel(typeLabel);
    classifResult->SetG(galaxyEvidence);
    classifResult->SetS(stellarEvidence);
    classifResult->SetQ(qsoEvidence);

    return;
}
