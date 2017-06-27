#ifndef lmfitcontroller_H
#define lmfitcontroller_H

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/ray.h>
#include <RedshiftLibrary/spectrum/template/template.h>

namespace NSEpic
{

class CLmfitController{
  public:
    CLmfitController(const CTemplate& tpl);
    bool isEmissionVelocityFitted();
    bool isAbsorptionVelocityFitted();
    bool isContinuumFitted();
    bool isLineTypeVelocityFitted(Int32 lineType);
    void resizeAmpsLine();
    void setAmpLine(Int32 idx, Float64 amp, Float64 ampErr);
    std::vector<Int32>  getFilteredIdx();
    bool removeNegAmpLine();
    void setNormFactor(Float64 normFactor);
    Float64 getNormFactor();
    bool needCalculateNormFactor();
    Int32 getNumberParameters();
    Int32 getIndEmissionVel();
    Int32 getIndAbsorptionVel();
    Int32 getIndContinuumAmp();
    Int32 getIndRedshift();
    void setContinummAmp(Float64 continuumAmp, Float64 continuumAmpErr);
    void setVelocityEmission(Float64 val, Float64 valErr);
    void setVelocityAbsorption(Float64 val, Float64 valErr);
    Float64 getAbsorptionVelocity();
    Float64 getEmissionVelocity();
    Float64 getContinuumAmp();
    void calculatedIndices();
    void addElement(Int32 elemId);
    void setMerit(Float64 merit);
    Float64 getMerit();
    Float64 getLineAmp(Int32 idx);
    Float64 getLineAmpErr(Int32 idx);
    bool isRedshiftFitted();
    CTemplate const * getTemplate();


  private:


    std::vector<Float64> m_ampsLinefitted;
    std::vector<Float64> m_ampErrLineFitted;
    bool m_normFactorSetted;
    Float64 m_normFactor;
    Float64 m_continuumAmp;
    Float64 m_continuumAmpErr;
    Float64 m_velAbs;
    Float64 m_velErrAbs;
    Float64 m_velEm;
    Float64 m_velErrEm;
    Float64 m_merit;
    std::vector<Int32> m_filteredEltsIdx;
    Int32 m_numberParameters;
    Int32 m_indAbsorptionVel;
    Int32 m_indEmissionVel;
    Int32 m_indContinuumAmp;
    Int32 m_indRedshift;
    bool m_NegAmpRemoved;
    CTemplate const * m_tpl;
};
}
#endif
