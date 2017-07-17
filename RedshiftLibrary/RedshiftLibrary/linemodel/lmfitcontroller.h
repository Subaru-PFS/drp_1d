#ifndef lmfitcontroller_H
#define lmfitcontroller_H

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/ray.h>
#include <RedshiftLibrary/spectrum/template/template.h>

namespace NSEpic
{

class CLmfitController{
  public:
    CLmfitController(
              const CTemplate& tpl,
              bool continumLoaded,
              bool continuumfit,
              bool emissionVelFit,
              bool absorptionVelFit
          );
    CLmfitController(
              bool emissionVelFit,
              bool absorptionVelFit
          );
    bool isEmissionVelocityFitted();
    bool isAbsorptionVelocityFitted();
    bool isContinuumFitted();
    bool isLineTypeVelocityFitted(Int32 lineType);
    bool isContinuumLoaded();
    bool isNoContinuum();
    void resizeAmpsLine();
    void setAmpLine(Int32 idx, Float64 amp, Float64 ampErr);
    std::vector<Int32>  getFilteredIdx();
    bool removeNegAmpLine();
    void setNormAmpLine(Float64 normAmpLine);
    Float64 getNormAmpLine();
    void setNormFactor(Float64 normFactor);
    Float64 getNormFactor();
    void setNormEmiFactor(Float64 normEmiFactor);
    Float64 getNormEmiFactor();
    void setNormAbsFactor(Float64 normAbsFactor);
    Float64 getNormAbsFactor();
    bool needCalculateNormFactor();
    Int32 getNumberParameters();
    Int32 getIndEmissionVel();
    Int32 getIndAbsorptionVel();
    Int32 getIndContinuumAmp();
    Int32 getIndRedshift();
    void setContinummAmp(Float64 continuumAmp, Float64 continuumAmpErr);
    void setVelocityEmission(Float64 val, Float64 valErr);
    void setVelocityAbsorption(Float64 val, Float64 valErr);
    void setRedshift(Float64 redshift, Float64 redshiftErr);
    Float64 getAbsorptionVelocity();
    Float64 getEmissionVelocity();
    Float64 getContinuumAmp();
    Float64 getRedshift();
    void calculatedIndices();
    void addElement(Int32 elemId);
    void setMerit(Float64 merit);
    Float64 getMerit();
    Float64 getLineAmp(Int32 idx);
    Float64 getLineAmpErr(Int32 idx);
    bool isRedshiftFitted();
    CTemplate const * getTemplate();


    Float64 lineAmp_LmToModel(Float64 lmLineAmp);
    Float64 lineAmp_ModelToLm(Float64 modelLineAmp);
    Float64 emiVel_LmToModel(Float64 lmEmiVel);
    Float64 emiVel_ModelToLm(Float64 modelEmiVel);
    Float64 absVel_LmToModel(Float64 lmAbsVel);
    Float64 absVel_ModelToLm(Float64 modelAbsVel);
    Float64 continuumAmp_LmToModel(Float64 lmContinuumAmp);
    Float64 continuumAmp_ModelToLm(Float64 lmContinuumAmp);

  private:

    bool m_continumLoaded ;
    bool m_continuumfit;
    bool m_noContinuum;
    bool m_emissionVelFit;
    bool m_absorptionVelFit ;
    std::vector<Float64> m_ampsLinefitted;
    std::vector<Float64> m_ampErrLineFitted;
    bool m_normFactorSetted;
    Float64 m_normFactor;
    Float64 m_normAmpLine;
    Float64 m_normEmiFactor;
    Float64 m_normAbsFactor;
    Float64 m_continuumAmp;
    Float64 m_continuumAmpErr;
    Float64 m_velAbs;
    Float64 m_velErrAbs;
    Float64 m_velEm;
    Float64 m_velErrEm;
    Float64 m_merit;
    Float64 m_redshift;
    Float64 m_redshiftErr;
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
