#include <RedshiftLibrary/linemodel/lmfitcontroller.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/ray.h>
#include <RedshiftLibrary/spectrum/template/template.h>

using namespace NSEpic;
/*
This class is use to drive  levenberg marquardt fitting.
Each instance correpond to one tentatie of fiiting. It has parameters of fitting, like which parameters need to be fit.
It's also store the results of the fitting, in order to temporize the moment when the model is set.
*/
CLmfitController::CLmfitController(
          const CTemplate& tpl,
          bool continumLoaded,
          bool continuumfit,
          bool emissionVelFit,
          bool absorptionVelFit
      ){
    m_tpl = &tpl;
    m_continumLoaded = continumLoaded;
    m_continuumfit = continuumfit;
    m_emissionVelFit = emissionVelFit;
    m_absorptionVelFit = absorptionVelFit;
}

CLmfitController::CLmfitController(
          bool emissionVelFit,
          bool absorptionVelFit
      ){
    m_tpl = NULL;
    m_continumLoaded = true;
    m_continuumfit = false;
    m_emissionVelFit = emissionVelFit;
    m_absorptionVelFit = absorptionVelFit;
}

bool CLmfitController::isEmissionVelocityFitted(){
  return m_emissionVelFit;
}

bool CLmfitController::isAbsorptionVelocityFitted(){
  return m_absorptionVelFit;
}

// is the continuum fitted by lmfit
bool CLmfitController::isContinuumFitted(){
  return m_continuumfit;
}

bool  CLmfitController::isContinuumLoaded(){
  return m_continumLoaded;
}

bool CLmfitController::isRedshiftFitted(){
  return true;;
}


bool CLmfitController::isLineTypeVelocityFitted(Int32 lineType){
  if(lineType ==CRay::nType_Emission){
    return isEmissionVelocityFitted();
  }else{
    return isAbsorptionVelocityFitted();
  }
}

void CLmfitController::resizeAmpsLine(){
  m_ampsLinefitted.resize(m_filteredEltsIdx.size());
  m_ampErrLineFitted.resize(m_filteredEltsIdx.size());
}


void CLmfitController::setAmpLine(Int32 idx, Float64 amp, Float64 ampErr){
  if(idx < m_ampsLinefitted.size()){
    m_ampsLinefitted[idx] = amp;
    m_ampErrLineFitted[idx] = ampErr;
  }
}

Float64  CLmfitController::getLineAmp(Int32 idx){
  if(idx < m_ampsLinefitted.size()){
    return m_ampsLinefitted[idx];
  }
  Log.LogError("idx greater than m_ampsLinefitted size");
  return -1.0;
}

Float64  CLmfitController::getLineAmpErr(Int32 idx){
  if(idx < m_ampErrLineFitted.size()){
    return m_ampErrLineFitted[idx];
  }
  Log.LogError("idx greater than m_ampErrLineFitted size");
  return -1.0;
}

void CLmfitController::addElement(Int32 elemId){
  m_filteredEltsIdx.push_back(elemId);
}

std::vector<Int32> CLmfitController::getFilteredIdx(){
  return m_filteredEltsIdx;
}

void CLmfitController::calculatedIndices(){
  m_numberParameters = m_filteredEltsIdx.size();

  if(isEmissionVelocityFitted()){
    m_indEmissionVel = m_numberParameters;
    m_numberParameters +=1;
  }
  if(isAbsorptionVelocityFitted()){
    m_indAbsorptionVel = m_numberParameters;
    m_numberParameters +=1;
  }
  if(isContinuumFitted()){
    m_indContinuumAmp = m_numberParameters;
    m_numberParameters +=1;
  }

  if(isRedshiftFitted()){
    m_indRedshift = m_numberParameters;
    m_numberParameters +=1;
  }

}

//======================================================================================
//Transform the variable from the model value to the value use in lmfit.
// if modif here modif the inverse function : modeltolm<-> lmtomodel and the df function in lmfitfunction (the jacobian is the derivided by lmvalue of the variable)
Float64 CLmfitController::lineAmp_LmToModel(Float64 lmLineAmp){
  return lmLineAmp*lmLineAmp;
  // return lmLineAmp;
}

Float64 CLmfitController::lineAmp_ModelToLm(Float64 modelLineAmp){
  return sqrt(modelLineAmp);
  // return modelLineAmp;
}

Float64 CLmfitController::emiVel_LmToModel(Float64 lmEmiVel){
  return lmEmiVel*lmEmiVel/m_normEmiFactor;
  // return lmEmiVel/m_normEmiFactor;
}

Float64 CLmfitController::emiVel_ModelToLm(Float64 modelEmiVel){
  return sqrt(modelEmiVel* m_normEmiFactor);
  // return modelEmiVel*m_normEmiFactor;
}

Float64 CLmfitController::absVel_LmToModel(Float64 lmAbsVel){
  return lmAbsVel*lmAbsVel/m_normAbsFactor;
  // return lmAbsVel/m_normAbsFactor;
}

Float64 CLmfitController::absVel_ModelToLm(Float64 modelAbsVel){
  return sqrt(modelAbsVel* m_normAbsFactor);
  // return modelAbsVel*m_normAbsFactor;
}

Float64 CLmfitController::continuumAmp_LmToModel(Float64 lmContinuumAmp){
  return lmContinuumAmp*lmContinuumAmp;
  // return lmContinuumAmp;
}

Float64 CLmfitController::continuumAmp_ModelToLm(Float64 lmContinuumAmp){
  return sqrt(lmContinuumAmp);
  // return lmContinuumAmp;
}

//======================================================================================
/*
Parse all the amplitude of element and look if some are neg. Those are remove from the fitting
return true if at least one line is removed
*/
bool CLmfitController::removeNegAmpLine(){
  m_NegAmpRemoved = false;
  for(Int32 ie=m_filteredEltsIdx.size()-1; ie>=0; ie--)
  {
      if(m_ampsLinefitted[ie]<0.0)
      {
          //Log.LogInfo( "LineModel Infos: erasing i= %d", ie);
          m_filteredEltsIdx.erase(m_filteredEltsIdx.begin() + ie);

          m_NegAmpRemoved = true;
      }
  }
  if(m_NegAmpRemoved){
    calculatedIndices();
  }
  return m_NegAmpRemoved;
}


void CLmfitController::setNormFactor(Float64 normFactor){
   m_normFactorSetted = true;
   m_normFactor = normFactor;
}
Float64 CLmfitController::getNormFactor(){
  return m_normFactor;
}

void CLmfitController::setNormEmiFactor(Float64 normEmiFactor){

   m_normEmiFactor = normEmiFactor;
}

Float64 CLmfitController::getNormEmiFactor(){
  return m_normEmiFactor;
}


void CLmfitController::setNormAbsFactor(Float64 normAbsFactor){

   m_normAbsFactor = normAbsFactor;
}

Float64 CLmfitController::getNormAbsFactor(){
  return m_normAbsFactor;
}


bool CLmfitController::needCalculateNormFactor(){
  if(!m_normFactorSetted){
    return true;
  }
  if(!isContinuumFitted() && m_NegAmpRemoved){
    return true;
  }
  return false;
}

Int32 CLmfitController::getNumberParameters(){
  return m_numberParameters;
}

Int32 CLmfitController::getIndEmissionVel(){
  return m_indEmissionVel;
}

Int32 CLmfitController::getIndAbsorptionVel(){
  return m_indAbsorptionVel;
}

Int32 CLmfitController::getIndContinuumAmp(){
  return m_indContinuumAmp;
}

Int32 CLmfitController::getIndRedshift(){
  return m_indRedshift;
}


void CLmfitController::setContinummAmp(Float64 continuumAmp, Float64 continuumAmpErr){
  m_continuumAmp = continuumAmp;
  m_continuumAmpErr = continuumAmpErr;
}

void  CLmfitController::setVelocityAbsorption(Float64 val, Float64 valErr){
  m_velAbs = val;
  m_velErrAbs = valErr;
}

void  CLmfitController::setVelocityEmission(Float64 val,Float64  valErr){
  m_velEm = val;
  m_velErrEm= valErr;
}

void CLmfitController::setRedshift(Float64 redshift, Float64 redshiftErr){
  m_redshift = redshift;
  m_redshiftErr = redshiftErr;
}

Float64  CLmfitController::getEmissionVelocity(){
  return m_velEm;
}

Float64  CLmfitController::getAbsorptionVelocity(){
  return m_velAbs;
}

Float64 CLmfitController::getContinuumAmp(){
  return m_continuumAmp;
}

Float64 CLmfitController::getRedshift(){
  return m_redshift;
}
void CLmfitController::setMerit(Float64 merit){
  m_merit = merit;
}
Float64 CLmfitController::getMerit(){
  return m_merit;
}

const CTemplate* CLmfitController::getTemplate(){
  return m_tpl;
}
