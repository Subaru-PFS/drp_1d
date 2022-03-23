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
#ifndef _REDSHIFT_LINEMODEL_LMFITCONTROLLER_
#define _REDSHIFT_LINEMODEL_LMFITCONTROLLER_

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/spectrum/template/template.h"

namespace NSEpic
{

class CLmfitController{
  public:
    CLmfitController(
              const std::shared_ptr<const CTemplate>& tpl,
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
    TInt32List  getFilteredIdx();
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
    Float64 getContinuumAmpErr();
    Float64 getRedshift();
    void calculatedIndices();
    void addElement(Int32 elemId);
    void setMerit(Float64 merit);
    Float64 getMerit();
    Float64 getLineAmp(Int32 idx);
    Float64 getLineAmpErr(Int32 idx);
    bool isRedshiftFitted();
    const std::shared_ptr<CTemplate const> getTemplate();


    Float64 lineAmp_LmToModel(Float64 lmLineAmp);
    Float64 lineAmp_ModelToLm(Float64 modelLineAmp);
    Float64 emiVel_LmToModel(Float64 lmEmiVel);
    Float64 emiVel_ModelToLm(Float64 modelEmiVel);
    Float64 absVel_LmToModel(Float64 lmAbsVel);
    Float64 absVel_ModelToLm(Float64 modelAbsVel);
    Float64 continuumAmp_LmToModel(Float64 lmContinuumAmp);
    Float64 continuumAmp_ModelToLm(Float64 lmContinuumAmp);
    //Float64 continuumAmpErr_LmToModel(Float64 lmContinuumAmpErr);
    //Float64 continuumAmpErr_ModelToLm(Float64 lmContinuumAmpErr);

  private:

    bool m_continumLoaded = false;
    bool m_continuumfit = false;
    bool m_noContinuum = false;
    bool m_emissionVelFit = false;
    bool m_absorptionVelFit = false;
    TFloat64List m_ampsLinefitted;
    TFloat64List m_ampErrLineFitted;
    bool m_normFactorSetted = false;
    Float64 m_normFactor = 0.;
    Float64 m_normAmpLine = 0.;
    Float64 m_normEmiFactor = 0.;
    Float64 m_normAbsFactor = 0.;
    Float64 m_continuumAmp = 0.;
    Float64 m_continuumAmpErr = 0.;
    Float64 m_velAbs = 0.;
    Float64 m_velErrAbs = 0.;
    Float64 m_velEm = 0.;
    Float64 m_velErrEm = 0.;
    Float64 m_merit = 0.;
    Float64 m_redshift = 0.;
    Float64 m_redshiftErr = 0.;
    TInt32List m_filteredEltsIdx;
    Int32 m_numberParameters = 0;
    Int32 m_indAbsorptionVel = 0;
    Int32 m_indEmissionVel = 0;
    Int32 m_indContinuumAmp = 0;
    Int32 m_indRedshift = 0;
    bool m_NegAmpRemoved = false;
    const std::shared_ptr<const CTemplate> m_tpl;
};
}
#endif
