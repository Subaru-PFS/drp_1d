#include <RedshiftLibrary/linemodel/lmfitcontroller.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/ray.h>
#include <RedshiftLibrary/spectrum/template/template.h>

using namespace NSEpic;
/*
This class is use to drive  levenberg marquardt fitting.
Each instance correpond to one tentatie of fiiting. It has parameters of fitting, like which parameters need to be fit.
It's also store the results of the fitting, in order to temporize the moment when the model is set. This allow to do several run of lmfit and
keep the best one.
*/
CLmfitController::CLmfitController(
          const CTemplate& tpl,
          bool continumLoaded,
          bool continuumfit,
          bool emissionVelFit,
          bool absorptionVelFit
      ){
    //m_tpl = &tpl;
    m_tpl = std::make_shared<const CTemplate>(tpl);
    m_noContinuum = false;
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
    m_noContinuum = true;
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

// return if the continuum is fitted by lmfit
bool CLmfitController::isContinuumFitted(){
  return m_continuumfit;
}

// is the continuum is loadded in linemodel ie interpolated on the grid.
bool  CLmfitController::isContinuumLoaded(){
  return m_continumLoaded;
}

// return if the line model don't not inclde a continuum
bool CLmfitController::isNoContinuum(){
    return m_noContinuum;
}

// return if the redshift is fitted by lmfit
bool CLmfitController::isRedshiftFitted(){
  return true;;
}

// return if the linetype (given in argument) is fitted by lmfit.
bool CLmfitController::isLineTypeVelocityFitted(Int32 lineType){
  if(lineType ==CRay::nType_Emission){
    return isEmissionVelocityFitted();
  }else{
    return isAbsorptionVelocityFitted();
  }
}

// resize the vectr of line amp and line error amp to the size of filtred element vector
void CLmfitController::resizeAmpsLine(){
  m_ampsLinefitted.resize(m_filteredEltsIdx.size());
  m_ampErrLineFitted.resize(m_filteredEltsIdx.size());
}

// store the value of amplitude of a line and it's error.
void CLmfitController::setAmpLine(Int32 idx, Float64 amp, Float64 ampErr){
  if(idx < m_ampsLinefitted.size()){
    m_ampsLinefitted[idx] = amp;
    m_ampErrLineFitted[idx] = ampErr;
  }
}

// return the stored amplitude of a element
Float64  CLmfitController::getLineAmp(Int32 idx){
  if(idx < m_ampsLinefitted.size()){
    return m_ampsLinefitted[idx];
  }
  Log.LogError("idx greater than m_ampsLinefitted size");
  return -1.0;
}

// return the stored error amplitude of a element
Float64  CLmfitController::getLineAmpErr(Int32 idx){
  if(idx < m_ampErrLineFitted.size()){
    return m_ampErrLineFitted[idx];
  }
  Log.LogError("idx greater than m_ampErrLineFitted size");
  return -1.0;
}

// add a element if the list of fiftted element
void CLmfitController::addElement(Int32 elemId){
  m_filteredEltsIdx.push_back(elemId);
}

// return the vecotr of index of fitted element
std::vector<UInt32> CLmfitController::getFilteredIdx(){
  return m_filteredEltsIdx;
}

// calculate the indice in the lmfit vector of varaible,
//it's take into account all the possibility of fit.(continuum or not, velocity or not ...)
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
  return lmLineAmp*lmLineAmp/m_normAmpLine;
  //return lmLineAmp;
}

Float64 CLmfitController::lineAmp_ModelToLm(Float64 modelLineAmp){
  return sqrt(modelLineAmp*m_normAmpLine);
  //return modelLineAmp;
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

/*Float64 CLmfitController::continuumAmpErr_LmToModel(Float64 lmContinuumAmpErr){
  return lmContinuumAmpErr*lmContinuumAmpErr;
  // return lmContinuumAmpErr;
}

Float64 CLmfitController::continuumAmpErr_ModelToLm(Float64 lmContinuumAmpErr){
  return sqrt(lmContinuumAmpErr);
  // return lmContinuumAmpErr;
}*/

//
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

// store the value of coef for element amplitude
void CLmfitController::setNormAmpLine(Float64 normAmpLineFactor){
    m_normAmpLine = normAmpLineFactor;
}

//return the value of coefficent factore for element amplitude
Float64 CLmfitController::getNormAmpLine(){
    return m_normAmpLine;
}

// store the coefficient that normalise the flux and the model.
void CLmfitController::setNormFactor(Float64 normFactor){
   m_normFactorSetted = true;
   m_normFactor = normFactor;
}

//return the coeffcient of normalisation for the flux an model
Float64 CLmfitController::getNormFactor(){
  return m_normFactor;
}

//store the coefficient of normalisation of the emission velocity factor
void CLmfitController::setNormEmiFactor(Float64 normEmiFactor){

   m_normEmiFactor = normEmiFactor;
}

//return the coefficient of normalisation of the emission velocity factor
Float64 CLmfitController::getNormEmiFactor(){
  return m_normEmiFactor;
}

//store the coefficient of normalisation of the absorption velocity factor
void CLmfitController::setNormAbsFactor(Float64 normAbsFactor){

   m_normAbsFactor = normAbsFactor;
}

//return the coefficient of normalisation of the absorption velocity factor
Float64 CLmfitController::getNormAbsFactor(){
  return m_normAbsFactor;
}

// return is the norma factor need to be recalculated.
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

// return the index of emission velocity in the vector of lmfit variable
Int32 CLmfitController::getIndEmissionVel(){
  return m_indEmissionVel;
}

// return the index of absorption velocity in the vector of lmfit variable
Int32 CLmfitController::getIndAbsorptionVel(){
  return m_indAbsorptionVel;
}

// return the index of continuum Amplitude in the vector of lmfit variable
Int32 CLmfitController::getIndContinuumAmp(){
  return m_indContinuumAmp;
}

// return the index of redshift in the vector of lmfit variable
Int32 CLmfitController::getIndRedshift(){
  return m_indRedshift;
}

//=========================================
// is this section with store and return the calculate value by lmfit.
// the value that is stored, is the value in linemodel.
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

Float64 CLmfitController::getContinuumAmpErr(){
  return m_continuumAmpErr;
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

std::shared_ptr<const CTemplate> CLmfitController::getTemplate(){
  return m_tpl;
}
