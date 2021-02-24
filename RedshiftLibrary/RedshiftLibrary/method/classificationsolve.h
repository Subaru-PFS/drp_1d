#ifndef _METHOD_CLASSIFICATION_
#define _METHOD_CLASSIFICATION_

#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/processflow/classificationresult.h>


namespace NSEpic
{

  class CClassificationSolve
  {

  public:

    CClassificationSolve(std::string enableStarFitting, std::string enableQsoFitting);
    ~CClassificationSolve();

    void Classify(std::shared_ptr<CSolveResult> galaxyResult, std::shared_ptr<CSolveResult> starResult, std::shared_ptr<CSolveResult> qspResult);
    
    std::string typeLabel = "U";//"G"/"S"/"Q"
    std::shared_ptr<CClassificationResult> classifResult;

  private:
  
    std::string m_enableStarFitting;
    std::string m_enableQsoFitting;
  };
}

#endif