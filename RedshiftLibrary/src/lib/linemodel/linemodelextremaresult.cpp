#include <RedshiftLibrary/linemodel/linemodelextremaresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <string>

using namespace NSEpic;


/**
 * \brief Empty constructor.
 **/
CLineModelExtremaResult::CLineModelExtremaResult()
{

}

/**
 * \brief Empty destructor.
 **/
CLineModelExtremaResult::~CLineModelExtremaResult()
{

}

void CLineModelExtremaResult::Resize(Int32 size)
{
    Extrema.resize(size);
    ExtremaMerit.resize(size);
    ExtremaMeritContinuum.resize(size);
    DeltaZ.resize(size);
    mTransposeM.resize(size);
    CorrScaleMarg.resize(size);
    NDof.resize(size);
    ExtremaLastPass.resize(size);
    lmfitPass.resize(size);
    snrHa.resize(size);
    lfHa.resize(size);
    snrOII.resize(size);
    lfOII.resize(size);

    Posterior.resize(size);
    StrongELSNR.resize(size);
    StrongELSNRAboveCut.resize(size);
    LogArea.resize(size);
    LogAreaCorrectedExtrema.resize(size);
    SigmaZ.resize(size);
    bic.resize(size);
    ContinuumIndexes.resize(size);
    OutsideLinesMask.resize(size);
    OutsideLinesSTDFlux.resize(size);
    OutsideLinesSTDError.resize(size);

    Elv.resize(size);
    Alv.resize(size);
    GroupsELv.resize(size);
    GroupsALv.resize(size);
    for(Int32 ke=0; ke<size; ke++)
    {
        GroupsELv[ke] = std::vector<Float64>(250, -1);   //WARNING: hardcoded ngroups max
        GroupsALv[ke] = std::vector<Float64>(250, -1);   //WARNING: hardcoded ngroups max
    }

    FittedTplName.resize(size);
    FittedTplAmplitude.resize(size);
    FittedTplMerit.resize(size);
    FittedTplDustCoeff.resize(size);
    FittedTplMeiksinIdx.resize(size);
    FittedTplRedshift.resize(size);
    FittedTplDtm.resize(size);
    FittedTplMtm.resize(size);
    FittedTplLogPrior.resize(size);
    FittedTplpCoeffs.resize(size);

    FittedTplshapeName.resize(size);
    FittedTplshapeIsmCoeff.resize(size);
    FittedTplshapeAmplitude.resize(size);
    FittedTplshapeDtm.resize(size);
    FittedTplshapeMtm.resize(size);
}

/**
 * \brief Empty method.
 **/
void CLineModelExtremaResult::SaveLine(const CDataStore &store, std::ostream& stream ) const
{

}

/**
 * \brief Prints the results currently in the argument store, in the argument stream.
 * Using the argument stream to print values from the argument store:
 * Print a header as a comment.
 * Print each redshift and chisquare values.
 * Print each Extrema as a comment.
 * Print each BIC as a comment.
 * Print each POSTERIOR as a comment.
 * Print each SigmaZ as a comment.
 * Print each LogArea as a comment.
 **/
void CLineModelExtremaResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    // save extrema list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#Extrema for z = {";
        for ( int i=0; i<Extrema.size(); i++)
        {
            stream <<  Extrema[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save extremaMerit list, on 1 line
    if(ExtremaMerit.size()>0){
        stream <<  "#ExtremaMerit for z = {";
        for ( int i=0; i<ExtremaMerit.size(); i++)
        {
            stream <<  ExtremaMerit[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save extremaMeritContinuum list, on 1 line
    if(ExtremaMeritContinuum.size()>0){
        stream <<  "#ExtremaMeritContinuum for z = {";
        for ( int i=0; i<ExtremaMeritContinuum.size(); i++)
        {
            stream <<  ExtremaMeritContinuum[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save extrema Deltaz list, on 1 line
    if(DeltaZ.size()>0){
        stream <<  "#ExtremaDeltaZ for z = {";
        for ( int i=0; i<DeltaZ.size(); i++)
        {
            stream <<  DeltaZ[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save extrema mTransposeM list, on 1 line
    if(mTransposeM.size()>0){
        stream <<  "#mTransposeM for z = {";
        for ( int i=0; i<mTransposeM.size(); i++)
        {
            stream <<  mTransposeM[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save NDof list, on 1 line
    if(NDof.size()>0){
        stream <<  "#NDof for z = {";
        for ( int i=0; i<NDof.size(); i++)
        {
            stream <<  NDof[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save CorrScaleMarg list, on 1 line
    if(CorrScaleMarg.size()>0){
        stream <<  "#CorrScaleMarg for z = {";
        for ( int i=0; i<CorrScaleMarg.size(); i++)
        {
            stream <<  CorrScaleMarg[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save ExtremaLastPass list, on 1 line
    if(ExtremaLastPass.size()>0){
        stream <<  "#ExtremaLastPass for z = {";
        for ( int i=0; i<ExtremaLastPass.size(); i++)
        {
            stream <<  ExtremaLastPass[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save lmfitPass list, on 1 line
    if(lmfitPass.size()>0){
        stream <<  "#lmfitPass for z = {";
        for ( int i=0; i<lmfitPass.size(); i++)
        {
            stream <<  lmfitPass[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save snrHa list, on 1 line
    if(snrHa.size()>0){
        stream <<  "#snrHa for z = {";
        for ( int i=0; i<snrHa.size(); i++)
        {
            stream <<  snrHa[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save lfHa list, on 1 line
    if(lfHa.size()>0){
        stream <<  "#lfHa for z = {";
        for ( int i=0; i<lfHa.size(); i++)
        {
            stream <<  lfHa[i] << "\t";
        }
        stream << "}" << std::endl;
    }


    // save snrOII list, on 1 line
    if(snrOII.size()>0){
        stream <<  "#snrOII for z = {";
        for ( int i=0; i<snrOII.size(); i++)
        {
            stream <<  snrOII[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save lfOII list, on 1 line
    if(lfOII.size()>0){
        stream <<  "#lfOII for z = {";
        for ( int i=0; i<lfOII.size(); i++)
        {
            stream <<  lfOII[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save bic list, on 1 line
    if(bic.size()>0){
        stream <<  "#BIC for each extrema = {";
        for ( int i=0; i<bic.size(); i++)
        {
            stream <<  bic[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save posterior list, on 1 line
    if(Posterior.size()>0){
        stream <<  "#POSTERIOR for each extrema = {";
        for ( int i=0; i<Posterior.size(); i++)
        {
            stream <<  Posterior[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save SigmaZ list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#SigmaZ for each extrema = {";
        for ( int i=0; i<SigmaZ.size(); i++)
        {
            stream <<  SigmaZ[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save LogArea list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#LogArea for each extrema = {";
        for ( int i=0; i<LogArea.size(); i++)
        {
            stream <<  LogArea[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save ContinuumIndexes list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#ContinuumIndexes Color for each extrema = {";
        for ( int i=0; i<ContinuumIndexes.size(); i++)
        {
            stream << "<";
            for(Int32 kci=0; kci<ContinuumIndexes[i].size(); kci++)
            {
                stream <<  ContinuumIndexes[i][kci].Color << "\t";
            }
            stream << ">";
        }
        stream << "}" << std::endl;
        stream <<  "#ContinuumIndexes Break for each extrema = {";
        for ( int i=0; i<ContinuumIndexes.size(); i++)
        {
            stream << "<";
            for(Int32 kci=0; kci<ContinuumIndexes[i].size(); kci++)
            {
                stream <<  ContinuumIndexes[i][kci].Break << "\t";
            }
            stream << ">";
        }
        stream << "}" << std::endl;
    }


    // save StrongELSNR list, on 1 line
    if(StrongELSNR.size()>0){
        stream <<  "#StrongELSNR for each extrema = {";
        for ( int i=0; i<StrongELSNR.size(); i++)
        {
            stream <<  StrongELSNR[i] << "\t";
        }
        stream << "}" << std::endl;
    }


    // save StrongELSNRAboveCut list, on 1 line
    if(StrongELSNRAboveCut.size()>0){
        stream <<  "#StrongELSNRAboveCut for each extrema = {";
        for ( int i=0; i<StrongELSNRAboveCut.size(); i++)
        {
            std::vector<std::string> line_list = StrongELSNRAboveCut[i];
            for ( int ki=0; ki<line_list.size(); ki++)
            {
                if(ki>0)
                {
                    stream << ",";
                }
                stream <<  line_list[ki];
            }
            stream << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplName, on 1 line
    if(FittedTplName.size()>0){
        stream <<  "#FittedTplName for each extrema = {";
        for ( int i=0; i<FittedTplName.size(); i++)
        {
            stream <<  FittedTplName[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplAmplitude, on 1 line
    if(FittedTplAmplitude.size()>0){
        stream <<  "#FittedTplAmplitude for each extrema = {";
        for ( int i=0; i<FittedTplAmplitude.size(); i++)
        {
            stream <<  FittedTplAmplitude[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplMerit, on 1 line
    if(FittedTplMerit.size()>0){
        stream <<  "#FittedTplMerit for each extrema = {";
        for ( int i=0; i<FittedTplMerit.size(); i++)
        {
            stream <<  FittedTplMerit[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplDustCoeff, on 1 line
    if(FittedTplDustCoeff.size()>0){
        stream <<  "#FittedTplDustCoeff for each extrema = {";
        stream << std::setprecision(3);
        for ( int i=0; i<FittedTplDustCoeff.size(); i++)
        {
            stream <<  FittedTplDustCoeff[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplMeiksinIdx, on 1 line
    if(FittedTplMeiksinIdx.size()>0){
        stream <<  "#FittedTplMeiksinIdx for each extrema = {";
        for ( int i=0; i<FittedTplMeiksinIdx.size(); i++)
        {
            stream <<  FittedTplMeiksinIdx[i] << "\t";
        }
        stream << "}" << std::endl;
    }


    // save FittedTplRedshift, on 1 line
    if(FittedTplRedshift.size()>0){
        stream <<  "#FittedTplRedshift for each extrema = {";
        stream << std::setprecision(8);
        for ( int i=0; i<FittedTplRedshift.size(); i++)
        {
            stream <<  FittedTplRedshift[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplDtm, on 1 line
    if(FittedTplDtm.size()>0){
        stream <<  "#FittedTplDtm for each extrema = {";
        stream << std::setprecision(8);
        for ( int i=0; i<FittedTplDtm.size(); i++)
        {
            stream <<  FittedTplDtm[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplMtm, on 1 line
    if(FittedTplMtm.size()>0){
        stream <<  "#FittedTplMtm for each extrema = {";
        stream << std::setprecision(8);
        for ( int i=0; i<FittedTplMtm.size(); i++)
        {
            stream <<  FittedTplMtm[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplLogPrior, on 1 line
    if(FittedTplLogPrior.size()>0){
        stream <<  "#FittedTplLogPrior for each extrema = {";
        stream << std::setprecision(8);
        for ( int i=0; i<FittedTplLogPrior.size(); i++)
        {
            stream <<  FittedTplLogPrior[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplpCoeffs, on 1 line
    if(FittedTplpCoeffs.size()>0){
        stream <<  "#FittedTplpCoeffs for each extrema = {";
        stream << std::setprecision(8);
        for ( int i=0; i<FittedTplpCoeffs.size(); i++)
        {
            if(i>0)
            {
                stream << "\t";
            }
            for ( int ip=0; ip<FittedTplpCoeffs[i].size(); ip++)
            {
                if(ip>0)
                {
                    stream << ":";
                }
                stream <<  FittedTplpCoeffs[i][ip];
            }

        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeName, on 1 line
    if(FittedTplshapeName.size()>0){
        stream <<  "#FittedTplshapeName for each extrema = {";
        for ( int i=0; i<FittedTplshapeName.size(); i++)
        {
            stream <<  FittedTplshapeName[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeIsmCoeff, on 1 line
    if(FittedTplshapeIsmCoeff.size()>0){
        stream <<  "#FittedTplshapeIsmCoeff for each extrema = {";
        for ( int i=0; i<FittedTplshapeIsmCoeff.size(); i++)
        {
            stream <<  FittedTplshapeIsmCoeff[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeAmplitude, on 1 line
    if(FittedTplshapeAmplitude.size()>0){
        stream <<  "#FittedTplshapeAmplitude for each extrema = {";
        for ( int i=0; i<FittedTplshapeAmplitude.size(); i++)
        {
            stream <<  FittedTplshapeAmplitude[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeDtm, on 1 line
    if(FittedTplshapeDtm.size()>0){
        stream <<  "#FittedTplshapeDtm for each extrema = {";
        for ( int i=0; i<FittedTplshapeDtm.size(); i++)
        {
            stream <<  FittedTplshapeDtm[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeMtm, on 1 line
    if(FittedTplshapeMtm.size()>0){
        stream <<  "#FittedTplshapeMtm for each extrema = {";
        for ( int i=0; i<FittedTplshapeMtm.size(); i++)
        {
            stream <<  FittedTplshapeMtm[i] << "\t";
        }
        stream << "}" << std::endl;
    }


    // save OutsideLinesSTDFlux, on 1 line
    if(OutsideLinesSTDFlux.size()>0){
        stream <<  "#OutsideLinesSTDFlux for each extrema = {";
        for ( int i=0; i<OutsideLinesSTDFlux.size(); i++)
        {
            stream << std::scientific << std::setprecision(5) <<  OutsideLinesSTDFlux[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save OutsideLinesSTDError, on 1 line
    if(OutsideLinesSTDError.size()>0){
        stream <<  "#OutsideLinesSTDError for each extrema = {";
        for ( int i=0; i<OutsideLinesSTDError.size(); i++)
        {
            stream << std::scientific << std::setprecision(5) <<  OutsideLinesSTDError[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save Elv, on 1 line
    if(Elv.size()>0){
        stream <<  "#Elv for each extrema = {";
        for ( int i=0; i<Elv.size(); i++)
        {
            stream << std::fixed << std::setprecision(1) <<  Elv[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save Alv, on 1 line
    if(Alv.size()>0){
        stream <<  "#Alv for each extrema = {";
        for ( int i=0; i<Alv.size(); i++)
        {
            stream << std::fixed << std::setprecision(1) <<  Alv[i] << "\t";
        }
        stream << "}" << std::endl;
    }



}

/**
 * \brief Prints the results currently in the argument store, in the argument stream.
 * Using the argument stream to print values from the argument store:
 * Print a header as a comment.
 * Print each redshift and chisquare values.
 * Print each Extrema as a comment.
 * Print each BIC as a comment.
 * Print each POSTERIOR as a comment.
 * Print each SigmaZ as a comment.
 * Print each LogArea as a comment.
 **/
void CLineModelExtremaResult::SaveJSON( const CDataStore& store, std::ostream& stream ) const
{
  stream << "{"<< std::endl;
  // save extrema list, on 1 line
  if(Extrema.size()>0){
    stream <<  "\"z_extrema\" : [";
    for ( int i=0; i<Extrema.size(); i++)
    {
      stream <<  Extrema[i];
      if( i< Extrema.size()-1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save extremaMerit list, on 1 line
  if(ExtremaMerit.size()>0){
    stream <<  "\"z_ExtremaMerit\" : [";
    for ( int i=0; i<ExtremaMerit.size(); i++)
    {
      stream << ExtremaMerit[i];
      if ( i<ExtremaMerit.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save extremaMeritContinuum list, on 1 line
  if(ExtremaMeritContinuum.size()>0){
    stream <<  "\"z_ExtremaMeritContinuum\" : [";
    for ( int i=0; i<ExtremaMeritContinuum.size(); i++)
    {
      stream << ExtremaMeritContinuum[i];
      if ( i<ExtremaMeritContinuum.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save extrema Deltaz list, on 1 line
  if(DeltaZ.size()>0){
    stream <<  "\"z_ExtremaDeltaZ\" : [";
    for ( int i=0; i<DeltaZ.size(); i++)
    {
      stream << DeltaZ[i];
      if ( i<DeltaZ.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save extrema mTransposeM list, on 1 line
  if(mTransposeM.size()>0){
    stream <<  "\"z_mTransposeM\" : [";
    for ( int i=0; i<mTransposeM.size(); i++)
    {
      stream << mTransposeM[i];
      if ( i<mTransposeM.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save NDof list, on 1 line
  if(NDof.size()>0){
    stream <<  "\"z_NDof\" : [";
    for ( int i=0; i<NDof.size(); i++)
    {
      stream << NDof[i];
      if ( i<NDof.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save CorrScaleMarg list, on 1 line
  if(CorrScaleMarg.size()>0){
    stream <<  "\"z_CorrScaleMarg\" : [";
    for ( int i=0; i<CorrScaleMarg.size(); i++)
    {
      stream << CorrScaleMarg[i];
      if ( i<CorrScaleMarg.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save ExtremaLastPass list, on 1 line
  if(ExtremaLastPass.size()>0){
    stream <<  "\"z_ExtremaLastPass\" : [";
    for ( int i=0; i<ExtremaLastPass.size(); i++)
    {
      stream << ExtremaLastPass[i];
      if ( i<ExtremaLastPass.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save lmfitPass list, on 1 line
  if(lmfitPass.size()>0){
    stream <<  "\"z_lmfitPass\" : [";
    for ( int i=0; i<lmfitPass.size(); i++)
    {
      stream << lmfitPass[i];
      if ( i<lmfitPass.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save snrHa list, on 1 line
  if(snrHa.size()>0){
    stream <<  "\"z_snrHa\" : [";
    for ( int i=0; i<snrHa.size(); i++)
    {
      stream << snrHa[i];
      if ( i<snrHa.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save lfHa list, on 1 line
  if(lfHa.size()>0){
    stream <<  "\"z_lfHa\" : [";
    for ( int i=0; i<lfHa.size(); i++)
    {
      stream << lfHa[i];
      if ( i<lfHa.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }


  // save snrOII list, on 1 line
  if(snrOII.size()>0){
    stream <<  "\"z_snrOII\" : [";
    for ( int i=0; i<snrOII.size(); i++)
    {
      stream << snrOII[i];
      if ( i<snrOII.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save lfOII list, on 1 line
  if(lfOII.size()>0){
    stream <<  "\"z_lfOII\" : [";
    for ( int i=0; i<lfOII.size(); i++)
    {
      stream << lfOII[i];
      if ( i<lfOII.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save bic list, on 1 line
  if(bic.size()>0){
    stream <<  "\"ext_BIC\" : [";
    for ( int i=0; i<bic.size(); i++)
    {
      stream << bic[i];
      if ( i<bic.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save posterior list, on 1 line
  if(Posterior.size()>0){
    stream <<  "\"ext_POSTERIOR\" : [";
    for ( int i=0; i<Posterior.size(); i++)
    {
      stream << Posterior[i];
      if ( i<Posterior.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save SigmaZ list, on 1 line
  if(Extrema.size()>0){
    stream <<  "\"ext_SigmaZ\" : [";
    for ( int i=0; i<SigmaZ.size(); i++)
    {
      stream << SigmaZ[i];
      if ( i<SigmaZ.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save LogArea list, on 1 line
  if(Extrema.size()>0){
    stream <<  "\"ext_LogArea\" : [";
    for ( int i=0; i<LogArea.size(); i++)
    {
      stream << LogArea[i];
      if ( i<LogArea.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save ContinuumIndexes list, on 1 line
  if(Extrema.size()>0){
    stream <<  "\"ext_ContinuumIndexesColor\" : [";
    for ( int i=0; i<ContinuumIndexes.size(); i++)
    {
      stream << "[";
      for(Int32 kci=0; kci<ContinuumIndexes[i].size(); kci++)
      {
	Float64 col = ContinuumIndexes[i][kci].Color;
        if (col != col ) stream << "null";
	else stream << col;
        if ( kci<ContinuumIndexes[i].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<ContinuumIndexes.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
    stream <<  "\"ext_ContinuumIndexesBreak\" :[";
    for ( int i=0; i<ContinuumIndexes.size(); i++)
    {
      stream << "[";
      for(Int32 kci=0; kci<ContinuumIndexes[i].size(); kci++)
      {
	Float64 break_ = ContinuumIndexes[i][kci].Break;
        if (break_ != break_ ) stream << "null";
	else stream <<  break_ ;
        if ( kci<ContinuumIndexes[i].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<ContinuumIndexes.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }


  // save StrongELSNR list, on 1 line
  if(StrongELSNR.size()>0){
    stream <<  "\"ext_StrongELSNR\" : [";
    for ( int i=0; i<StrongELSNR.size(); i++)
    {
      stream <<StrongELSNR[i];
      if ( i<StrongELSNR.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }


  // save StrongELSNRAboveCut list, on 1 line
  if(StrongELSNRAboveCut.size()>0){
    stream <<  "\"ext_StrongELSNRAboveCut\" : [";
    for ( int i=0; i<StrongELSNRAboveCut.size(); i++)
    {
      std::vector<std::string> line_list = StrongELSNRAboveCut[i];
      stream << "[";
      for ( int ki=0; ki<line_list.size(); ki++)
      {
        stream <<  "\"" << line_list[ki] << "\"";
        if (ki < line_list.size()-1) stream <<",";
      }
      stream << "]";
      if (i<StrongELSNRAboveCut.size()-1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplName, on 1 line
  if(FittedTplName.size()>0){
    stream <<  "\"ext_FittedTplName\" : [";
    for ( int i=0; i<FittedTplName.size(); i++)
    {
      stream << "\""<<FittedTplName[i]<<"\"";
      if ( i<FittedTplName.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplAmplitude, on 1 line
  if(FittedTplAmplitude.size()>0){
    stream <<  "\"ext_FittedTplAmplitude\" : [";
    for ( int i=0; i<FittedTplAmplitude.size(); i++)
    {
      stream << FittedTplAmplitude[i];
      if ( i<FittedTplAmplitude.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplMerit, on 1 line
  if(FittedTplMerit.size()>0){
    stream <<  "\"ext_FittedTplMerit\" : [";
    for ( int i=0; i<FittedTplMerit.size(); i++)
    {
      stream << FittedTplMerit[i];
      if ( i<FittedTplMerit.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplDustCoeff, on 1 line
  if(FittedTplDustCoeff.size()>0){
    stream <<  "\"ext_FittedTplDustCoeff\" : [";
    stream << std::setprecision(3);
    for ( int i=0; i<FittedTplDustCoeff.size(); i++)
    {
      stream << FittedTplDustCoeff[i];
      if ( i<FittedTplDustCoeff.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplMeiksinIdx, on 1 line
  if(FittedTplMeiksinIdx.size()>0){
    stream <<  "\"ext_FittedTplMeiksinIdx\" : [";
    for ( int i=0; i<FittedTplMeiksinIdx.size(); i++)
    {
      stream << FittedTplMeiksinIdx[i];
      if ( i<FittedTplMeiksinIdx.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }


  // save FittedTplRedshift, on 1 line
  if(FittedTplRedshift.size()>0){
    stream <<  "\"ext_FittedTplRedshift\" : [";
    stream << std::setprecision(8);
    for ( int i=0; i<FittedTplRedshift.size(); i++)
    {
      stream << FittedTplRedshift[i];
      if ( i<FittedTplRedshift.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplDtm, on 1 line
  if(FittedTplDtm.size()>0){
    stream <<  "\"ext_FittedTplDtm\" : [";
    stream << std::setprecision(8);
    for ( int i=0; i<FittedTplDtm.size(); i++)
    {
      stream << FittedTplDtm[i];
      if ( i<FittedTplDtm.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplMtm, on 1 line
  if(FittedTplMtm.size()>0){
    stream <<  "\"ext_FittedTplMtm\" : [";
    stream << std::setprecision(8);
    for ( int i=0; i<FittedTplMtm.size(); i++)
    {
      stream << FittedTplMtm[i];
      if ( i<FittedTplMtm.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplLogPrior, on 1 line
  if(FittedTplLogPrior.size()>0){
    stream <<  "\"ext_FittedTplLogPrior\" : [";
    stream << std::setprecision(8);
    for ( int i=0; i<FittedTplLogPrior.size(); i++)
    {
      stream << FittedTplLogPrior[i];
      if ( i<FittedTplLogPrior.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplpCoeffs, on 1 line
  if(FittedTplpCoeffs.size()>0){
    stream <<  "\"ext_FittedTplpCoeffs\" : [";
    stream << std::setprecision(8);
    for ( int i=0; i<FittedTplpCoeffs.size(); i++)
    {
      stream << "[";
      for ( int ip=0; ip<FittedTplpCoeffs[i].size(); ip++)
      {
        stream <<  FittedTplpCoeffs[i][ip];
        if ( ip< FittedTplpCoeffs[i].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<FittedTplpCoeffs.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplshapeName, on 1 line
  if(FittedTplshapeName.size()>0){
    stream <<  "\"ext_FittedTplshapeName\" : [";
    for ( int i=0; i<FittedTplshapeName.size(); i++)
    {
      stream << "\"" << FittedTplshapeName[i] << "\"";
      if ( i<FittedTplshapeName.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplshapeIsmCoeff, on 1 line
  if(FittedTplshapeIsmCoeff.size()>0){
    stream <<  "\"ext_FittedTplshapeIsmCoeff\" : [";
    for ( int i=0; i<FittedTplshapeIsmCoeff.size(); i++)
    {
      stream << FittedTplshapeIsmCoeff[i];
      if ( i<FittedTplshapeIsmCoeff.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplshapeAmplitude, on 1 line
  if(FittedTplshapeAmplitude.size()>0){
    stream <<  "\"ext_FittedTplshapeAmplitude\" : [";
    for ( int i=0; i<FittedTplshapeAmplitude.size(); i++)
    {
      stream << FittedTplshapeAmplitude[i];
      if ( i<FittedTplshapeAmplitude.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplshapeDtm, on 1 line
  if(FittedTplshapeDtm.size()>0){
    stream <<  "\"ext_FittedTplshapeDtm\" : [";
    for ( int i=0; i<FittedTplshapeDtm.size(); i++)
    {
      stream << FittedTplshapeDtm[i];
      if ( i<FittedTplshapeDtm.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save FittedTplshapeMtm, on 1 line
  if(FittedTplshapeMtm.size()>0){
    stream <<  "\"ext_FittedTplshapeMtm\" : [";
    for ( int i=0; i<FittedTplshapeMtm.size(); i++)
    {
      stream << FittedTplshapeMtm[i];
      if ( i<FittedTplshapeMtm.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }


  // save OutsideLinesSTDFlux, on 1 line
  if(OutsideLinesSTDFlux.size()>0){
    stream <<  "\"ext_OutsideLinesSTDFlux\" : [";
    for ( int i=0; i<OutsideLinesSTDFlux.size(); i++)
    {
      stream << std::scientific << std::setprecision(5) <<  OutsideLinesSTDFlux[i];
      if ( i<OutsideLinesSTDFlux.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save OutsideLinesSTDError, on 1 line
  if(OutsideLinesSTDError.size()>0){
    stream <<  "\"ext_OutsideLinesSTDError\" : [";
    for ( int i=0; i<OutsideLinesSTDError.size(); i++)
    {
      stream << std::scientific << std::setprecision(5) <<  OutsideLinesSTDError[i] ;
      if ( i<OutsideLinesSTDError.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save Elv, on 1 line
  if(Elv.size()>0){
    stream <<  "\"ext_Elv\" : [";
    for ( int i=0; i<Elv.size(); i++)
    {
      stream << std::fixed << std::setprecision(1) <<  Elv[i];
      if ( i<Elv.size() - 1) stream << ",";
    }
    stream << "]," << std::endl;
  }

  // save Alv, on 1 line
  if(Alv.size()>0){
    stream <<  "\"ext_Alv\" : [";
    for ( int i=0; i<Alv.size(); i++)
    {
      stream << std::fixed << std::setprecision(1) <<  Alv[i] ;
      if ( i<Alv.size() - 1) stream << ",";
    }
    stream << "]"<< std::endl << "}" << std::endl;
  }



}

