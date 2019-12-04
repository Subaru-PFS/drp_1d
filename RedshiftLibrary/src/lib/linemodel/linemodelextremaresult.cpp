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
  SaveTFloat64List(stream,"z_extrema", Extrema);
  stream << "," << std::endl;

  // save extremaMerit list, on 1 line
  SaveTFloat64List(stream,"z_ExtremaMerit",ExtremaMerit);
  stream << "," << std::endl;
  // save extremaMeritContinuum list, on 1 line
  SaveTFloat64List(stream,"z_ExtremaMeritContinuum",ExtremaMeritContinuum);
  stream << "," << std::endl;
  // save extrema Deltaz list, on 1 line
  SaveTFloat64List(stream,"z_ExtremaDeltaZ",DeltaZ);
  stream << "," << std::endl;
  // save extrema mTransposeM list, on 1 line
  SaveTFloat64List(stream,"z_mTransposeM",mTransposeM);
  stream << "," << std::endl;
  // save NDof list, on 1 line
  SaveInt32Vector(stream,"z_NDof",NDof);
  stream << "," << std::endl;
  // save CorrScaleMarg list, on 1 line
  SaveTFloat64List(stream,"z_CorrScaleMarg",CorrScaleMarg);
  stream << "," << std::endl;
  // save ExtremaLastPass list, on 1 line
  SaveTFloat64List(stream,"z_ExtremaLastPass",ExtremaLastPass);
  stream << "," << std::endl;
  // save lmfitPass list, on 1 line
  SaveTFloat64List(stream,"z_lmfitPass",lmfitPass);
  stream << "," << std::endl;
  // save snrHa list, on 1 line
  SaveTFloat64List(stream,"z_snrHa",snrHa);
  stream << "," << std::endl;
  // save lfHa list, on 1 line
  SaveTFloat64List(stream,"z_lfHa",lfHa);
  stream << "," << std::endl;

  // save snrOII list, on 1 line
  SaveTFloat64List(stream,"z_snrOII",snrOII);
  stream << "," << std::endl;
  // save lfOII list, on 1 line
  SaveTFloat64List(stream,"z_lfOII",lfOII);
  stream << "," << std::endl;
  // save bic list, on 1 line
  SaveTFloat64List(stream,"ext_BIC",bic);
  stream << "," << std::endl;
  // save posterior list, on 1 line
  SaveTFloat64List(stream,"ext_POSTERIOR",Posterior);
  stream << "," << std::endl;
  // save SigmaZ list, on 1 line
  SaveTFloat64List(stream,"ext_SigmaZ",SigmaZ);
  stream << "," << std::endl;
  // save LogArea list, on 1 line
  SaveTFloat64List(stream,"ext_LogArea",LogArea);
  stream << "," << std::endl;
  // save ContinuumIndexes list, on 1 line
  SaveTContinuumIndexListVector(stream,"ext_ContinuumIndexes",ContinuumIndexes);
  stream << "," << std::endl;	


  // save StrongELSNR list, on 1 line
  SaveTFloat64List(stream,"ext_StrongELSNR",StrongELSNR);
  stream << "," << std::endl;

  // save StrongELSNRAboveCut list, on 1 line
  SaveStringVectorOfVector(stream,"ext_StrongELSNRAboveCut",StrongELSNRAboveCut);
  stream << "," << std::endl;

  // save FittedTplName, on 1 line
  SaveStringVector(stream,"ext_FittedTplName",FittedTplName);
  stream << "," << std::endl;
  // save FittedTplAmplitude, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplAmplitude",FittedTplAmplitude);
  stream << "," << std::endl;
  // save FittedTplMerit, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplMerit",FittedTplMerit);
  stream << "," << std::endl;
  // save FittedTplDustCoeff, on 1 line
  stream << std::setprecision(3);
  SaveTFloat64List(stream,"ext_FittedTplDustCoeff",FittedTplDustCoeff); 
  stream << "," << std::endl; 

  // save FittedTplMeiksinIdx, on 1 line
  SaveInt32Vector(stream,"ext_FittedTplMeiksinIdx",FittedTplMeiksinIdx);
  stream << "," << std::endl;
   
  // save FittedTplRedshift, on 1 line
   stream << std::setprecision(8);
  SaveTFloat64List(stream,"ext_FittedTplRedshift",FittedTplRedshift);  
  stream << "," << std::endl;

  // save FittedTplDtm, on 1 line
      stream << std::setprecision(8);
  SaveTFloat64List(stream,"ext_FittedTplDtm",FittedTplDtm);
  stream << "," << std::endl;
  // save FittedTplMtm, on 1 line
    stream << std::setprecision(8);
  SaveTFloat64List(stream,"ext_FittedTplMtm",FittedTplMtm);
  stream << "," << std::endl;
  // save FittedTplLogPrior, on 1 line
    stream << std::setprecision(8);
  SaveTFloat64List(stream,"ext_FittedTplLogPrior",FittedTplLogPrior);
  stream << "," << std::endl;
  // save FittedTplpCoeffs, on 1 line
    stream << std::setprecision(8);
  SaveTFloat64ListOfList(stream,"ext_FittedTplpCoeffs",FittedTplpCoeffs);
  stream << "," << std::endl;

  // save FittedTplshapeName, on 1 line
  SaveStringVector(stream,"ext_FittedTplshapeName",FittedTplshapeName);
  stream << "," << std::endl;
  // save FittedTplshapeIsmCoeff, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplshapeIsmCoeff",FittedTplshapeIsmCoeff);
  stream << "," << std::endl;
  // save FittedTplshapeAmplitude, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplshapeAmplitude",FittedTplshapeAmplitude);
  stream << "," << std::endl;
  // save FittedTplshapeDtm, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplshapeDtm",FittedTplshapeDtm);
  stream << "," << std::endl;
  // save FittedTplshapeMtm, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplshapeMtm",FittedTplshapeMtm);
  stream << "," << std::endl;

  // save OutsideLinesSTDFlux, on 1 line
  stream << std::scientific << std::setprecision(5);
  SaveTFloat64List(stream,"ext_OutsideLinesSTDFlux",OutsideLinesSTDFlux);
  stream << "," << std::endl;
  // save OutsideLinesSTDError, on 1 line
  stream << std::scientific << std::setprecision(5);
  SaveTFloat64List(stream,"ext_OutsideLinesSTDError",OutsideLinesSTDError);
  stream << "," << std::endl;
  // save Elv, on 1 line
  stream << std::fixed << std::setprecision(1);
  SaveTFloat64List(stream,"ext_Elv",Elv);
  stream << "," << std::endl;
  // save Alv, on 1 line
  stream << std::fixed << std::setprecision(1);
  SaveTFloat64List(stream,"ext_Alv",Alv);

  stream <<  "}";

}

