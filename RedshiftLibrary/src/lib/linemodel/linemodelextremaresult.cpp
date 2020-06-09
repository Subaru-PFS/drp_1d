#include <RedshiftLibrary/linemodel/linemodelextremaresult.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/log/log.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <iostream>
#include <numeric>
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
    ExtremaPDF.resize(size);
    ExtremaIDs.resize(size);
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
    FittedTplAmplitudeError.resize(size);
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

<<<<<<< HEAD
<<<<<<< HEAD

bool CLineModelExtremaResult::RemoveSecondPassCandidatebyIdx(Int32 idx){
    ExtremaPDF.erase(ExtremaPDF.begin() + idx);
    ExtremaIDs.erase( ExtremaIDs.begin() + idx);
    Extrema.erase( Extrema.begin() + idx);
    ExtremaMerit.erase( ExtremaMerit.begin() + idx);
    ExtremaMeritContinuum.erase( ExtremaMeritContinuum.begin() + idx);
    DeltaZ.erase( DeltaZ.begin() + idx);
    mTransposeM.erase( mTransposeM.begin() + idx);
    CorrScaleMarg.erase( CorrScaleMarg.begin() + idx);
    NDof.erase( NDof.begin() + idx);
    ExtremaLastPass.erase( ExtremaLastPass.begin() + idx);
    lmfitPass.erase( lmfitPass.begin() + idx);
    snrHa.erase( snrHa.begin() + idx);
    lfHa.erase( lfHa.begin() + idx);
    snrOII.erase( snrOII.begin() + idx);
    lfOII.erase( lfOII.begin() + idx);

    Posterior.erase( Posterior.begin() + idx);
    StrongELSNR.erase( StrongELSNR.begin() + idx);
    StrongELSNRAboveCut.erase( StrongELSNRAboveCut.begin() + idx);
    LogArea.erase( LogArea.begin() + idx);
    LogAreaCorrectedExtrema.erase( LogAreaCorrectedExtrema.begin() + idx);
    SigmaZ.erase( SigmaZ.begin() + idx);
    bic.erase( bic.begin() + idx);
    ContinuumIndexes.erase( ContinuumIndexes.begin() + idx);
    OutsideLinesMask.erase( OutsideLinesMask.begin() + idx);
    OutsideLinesSTDFlux.erase( OutsideLinesSTDFlux.begin() + idx);
    OutsideLinesSTDError.erase( OutsideLinesSTDError.begin() + idx);

    Elv.erase( Elv.begin() + idx);
    Alv.erase( Alv.begin() + idx);
    GroupsELv.erase( GroupsELv.begin() + idx);
    GroupsALv.erase( GroupsALv.begin() + idx);

    FittedTplName.erase( FittedTplName.begin() + idx);
    FittedTplAmplitude.erase( FittedTplAmplitude.begin() + idx);
    FittedTplMerit.erase( FittedTplMerit.begin() + idx);
    FittedTplDustCoeff.erase( FittedTplDustCoeff.begin() + idx);
    FittedTplMeiksinIdx.erase( FittedTplMeiksinIdx.begin() + idx);
    FittedTplRedshift.erase( FittedTplRedshift.begin() + idx);
    FittedTplDtm.erase( FittedTplDtm.begin() + idx);
    FittedTplMtm.erase( FittedTplMtm.begin() + idx);
    FittedTplLogPrior.erase( FittedTplLogPrior.begin() + idx);
    FittedTplpCoeffs.erase( FittedTplpCoeffs.begin() + idx);

    FittedTplshapeName.erase( FittedTplshapeName.begin() + idx);
    FittedTplshapeIsmCoeff.erase( FittedTplshapeIsmCoeff.begin() + idx);
    FittedTplshapeAmplitude.erase( FittedTplshapeAmplitude.begin() + idx);
    FittedTplshapeDtm.erase( FittedTplshapeDtm.begin() + idx);
    FittedTplshapeMtm.erase(FittedTplshapeMtm.begin() + idx);
    ExtremaExtendedRedshifts.erase(ExtremaExtendedRedshifts.begin() + idx);
    return true;
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
void CLineModelExtremaResult::Save( const CDataStore& store, std::ostream& stream) const
{
    //here there is no scope set--> error in accessing results using GetScope()
    auto res = store.GetGlobalResult( "candidatesresult"  );
    auto candResults = std::dynamic_pointer_cast<const CPdfCandidateszResult>( res.lock());
    TInt32List Rank_PDF = candResults->Rank;

    //below is a sign that we are saving firstpass data
    //TODO: check if removing these info is valid for all types of spectra
    bool zeros = std::all_of(ExtremaMeritContinuum.begin(), ExtremaMeritContinuum.end(), [](int i) { return i==0; });    
    if(zeros){
         Int32 n = Rank_PDF.size();
        if(n<Extrema.size()){
            //in this case, rank doesnt mean anything
            /*Rank_PDF.clear();
            for(Int32 i = 0; i<Extrema.size(); i++){
                Rank_PDF.push_back(i);
            }*/
            //in this case, rank doesnt mean anything
            Rank_PDF.clear();
            for(Int32 i = 0; i<Extrema.size(); i++){
                Rank_PDF.push_back(i);
            }
            std::vector<std::string> ids = candResults->ExtremaIDs;
            Int32 ret = FixRanksUsingSortedIDs(Rank_PDF, ids);
            if(ret<0){
                Log.LogError("   CLineModelExtremaResult::FixRanksUsingSortedIDs failed!");
                throw std::runtime_error("   CLineModelExtremaResult::FixRanksUsingSortedIDs failed!");
            }

        }
    }

    // save extrema list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#Extrema for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  Extrema[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }
    if(ExtremaIDs.size()>0){
        stream <<  "#ExtremaIDs for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            std::string s = ExtremaIDs[Rank_PDF[i]];
            stream <<  s  << "\t";
        }
        stream << "}" << std::endl;
    }

    if(!zeros){
        // save extrema reference rank, used to map
        if(Extrema.size()>0){
          stream <<  "#Extrema final Rank for z = {";
            for ( int i=0; i<Rank_PDF.size(); i++)
            {
                stream <<  i << "\t";
            }
            stream << "}" << std::endl;
        }

     TFloat64List ExtremaPDF = candResults->ValSumProba;
     if(ExtremaPDF.size()>0){
        stream <<  "#Extrema IntgPDF for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  ExtremaPDF[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
     }
    }
    // save extremaMerit list, on 1 line
    if(ExtremaMerit.size()>0){
        stream <<  "#ExtremaMerit for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  ExtremaMerit[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }
if(!zeros){
    // save extremaMeritContinuum list, on 1 line
    if(ExtremaMeritContinuum.size()>0){
        stream <<  "#ExtremaMeritContinuum for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  ExtremaMeritContinuum[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save extrema Deltaz list, on 1 line
    if(DeltaZ.size()>0){
        stream <<  "#ExtremaDeltaZ for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  DeltaZ[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save extrema mTransposeM list, on 1 line
    if(mTransposeM.size()>0){
        stream <<  "#mTransposeM for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  mTransposeM[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save NDof list, on 1 line
    if(NDof.size()>0){
        stream <<  "#NDof for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  NDof[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save CorrScaleMarg list, on 1 line
    if(CorrScaleMarg.size()>0){
        stream <<  "#CorrScaleMarg for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  CorrScaleMarg[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save ExtremaLastPass list, on 1 line
    if(ExtremaLastPass.size()>0){
        stream <<  "#ExtremaLastPass for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  ExtremaLastPass[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save lmfitPass list, on 1 line
    if(lmfitPass.size()>0){
        stream <<  "#lmfitPass for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  lmfitPass[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save snrHa list, on 1 line
    if(snrHa.size()>0){
        stream <<  "#snrHa for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  snrHa[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save lfHa list, on 1 line
    if(lfHa.size()>0){
        stream <<  "#lfHa for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  lfHa[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }


    // save snrOII list, on 1 line
    if(snrOII.size()>0){
        stream <<  "#snrOII for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  snrOII[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save lfOII list, on 1 line
    if(lfOII.size()>0){
        stream <<  "#lfOII for z = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  lfOII[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save bic list, on 1 line
    if(bic.size()>0){
        stream <<  "#BIC for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  bic[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save posterior list, on 1 line
    if(Posterior.size()>0){
        stream <<  "#POSTERIOR for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  Posterior[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save SigmaZ list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#SigmaZ for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  SigmaZ[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save LogArea list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#LogArea for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  LogArea[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save ContinuumIndexes list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#ContinuumIndexes Color for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream << "<";
            for(Int32 kci=0; kci<ContinuumIndexes[Rank_PDF[i]].size(); kci++)
            {
                stream <<  ContinuumIndexes[Rank_PDF[i]][kci].Color << "\t";
            }
            stream << ">";
        }
        stream << "}" << std::endl;
        stream <<  "#ContinuumIndexes Break for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream << "<";
            for(Int32 kci=0; kci<ContinuumIndexes[Rank_PDF[i]].size(); kci++)
            {
                stream <<  ContinuumIndexes[Rank_PDF[i]][kci].Break << "\t";
            }
            stream << ">";
        }
        stream << "}" << std::endl;
    }


    // save StrongELSNR list, on 1 line
    if(StrongELSNR.size()>0){
        stream <<  "#StrongELSNR for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  StrongELSNR[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }


    // save StrongELSNRAboveCut list, on 1 line
    if(StrongELSNRAboveCut.size()>0){
        stream <<  "#StrongELSNRAboveCut for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            std::vector<std::string> line_list = StrongELSNRAboveCut[Rank_PDF[i]];
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
}
    // save FittedTplName, on 1 line
    if(FittedTplName.size()>0){
        stream <<  "#FittedTplName for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplName[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplAmplitude, on 1 line
    if(FittedTplAmplitude.size()>0){
        stream <<  "#FittedTplAmplitude for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplAmplitude[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplAmplitudeError, on 1 line
    if(FittedTplAmplitudeError.size()>0){
        stream <<  "#FittedTplAmplitudeError for each extrema = {";
        for ( int i=0; i<FittedTplAmplitudeError.size(); i++)
        {
            stream <<  FittedTplAmplitudeError[Rank_PDF[i]] << "\t";
        }
	stream << "}" << std::endl;
    }

    // save FittedTplMerit, on 1 line
    if(FittedTplMerit.size()>0){
        stream <<  "#FittedTplMerit for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplMerit[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplDustCoeff, on 1 line
    if(FittedTplDustCoeff.size()>0){
        stream <<  "#FittedTplDustCoeff for each extrema = {";
        stream << std::setprecision(3);
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplDustCoeff[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplMeiksinIdx, on 1 line
    if(FittedTplMeiksinIdx.size()>0){
        stream <<  "#FittedTplMeiksinIdx for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplMeiksinIdx[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }


    // save FittedTplRedshift, on 1 line
    if(FittedTplRedshift.size()>0){
        stream <<  "#FittedTplRedshift for each extrema = {";
        stream << std::setprecision(8);
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplRedshift[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplDtm, on 1 line
    if(FittedTplDtm.size()>0){
        stream <<  "#FittedTplDtm for each extrema = {";
        stream << std::setprecision(8);
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplDtm[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplMtm, on 1 line
    if(FittedTplMtm.size()>0){
        stream <<  "#FittedTplMtm for each extrema = {";
        stream << std::setprecision(8);
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplMtm[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }
 if(!zeros){
    // save FittedTplLogPrior, on 1 line
    if(FittedTplLogPrior.size()>0){
        stream <<  "#FittedTplLogPrior for each extrema = {";
        stream << std::setprecision(8);
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplLogPrior[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplpCoeffs, on 1 line
    if(FittedTplpCoeffs.size()>0){
        stream <<  "#FittedTplpCoeffs for each extrema = {";
        stream << std::setprecision(8);
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            if(i>0)
            {
                stream << "\t";
            }
            for ( int ip=0; ip<FittedTplpCoeffs[Rank_PDF[i]].size(); ip++)
            {
                if(ip>0)
                {
                    stream << ":";
                }
                stream <<  FittedTplpCoeffs[Rank_PDF[i]][ip];
            }

        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeName, on 1 line
    if(FittedTplshapeName.size()>0){
        stream <<  "#FittedTplshapeName for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplshapeName[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeIsmCoeff, on 1 line
    if(FittedTplshapeIsmCoeff.size()>0){
        stream <<  "#FittedTplshapeIsmCoeff for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplshapeIsmCoeff[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeAmplitude, on 1 line
    if(FittedTplshapeAmplitude.size()>0){
        stream <<  "#FittedTplshapeAmplitude for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplshapeAmplitude[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeDtm, on 1 line
    if(FittedTplshapeDtm.size()>0){
        stream <<  "#FittedTplshapeDtm for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplshapeDtm[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeMtm, on 1 line
    if(FittedTplshapeMtm.size()>0){
        stream <<  "#FittedTplshapeMtm for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream <<  FittedTplshapeMtm[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }


    // save OutsideLinesSTDFlux, on 1 line
    if(OutsideLinesSTDFlux.size()>0){
        stream <<  "#OutsideLinesSTDFlux for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream << std::scientific << std::setprecision(5) <<  OutsideLinesSTDFlux[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save OutsideLinesSTDError, on 1 line
    if(OutsideLinesSTDError.size()>0){
        stream <<  "#OutsideLinesSTDError for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream << std::scientific << std::setprecision(5) <<  OutsideLinesSTDError[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }
 }
    // save Elv, on 1 line
    if(Elv.size()>0){
        stream <<  "#Elv for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream << std::fixed << std::setprecision(1) <<  Elv[Rank_PDF[i]] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save Alv, on 1 line
    if(Alv.size()>0){
        stream <<  "#Alv for each extrema = {";
        for ( int i=0; i<Rank_PDF.size(); i++)
        {
            stream << std::fixed << std::setprecision(1) <<  Alv[Rank_PDF[i]] << "\t";
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
void CLineModelExtremaResult::SaveJSON( const CDataStore& store, std::ostream& stream) const
{
  auto res = store.GetGlobalResult( "candidatesresult");
  auto candResults = std::dynamic_pointer_cast<const CPdfCandidateszResult>( res.lock());
  TInt32List order = candResults->Rank;

<<<<<<< HEAD
  bool zeros = std::all_of(ExtremaMeritContinuum.begin(), ExtremaMeritContinuum.end(), [](int i) { return i==0; });
  if(zeros){
        //saving firstpass data
        //in the case where some candidates from first-pass are eliminated from second-pass
        //make sure to save all data we collected from first pass
        Int32 n = order.size();
        if(n<Extrema.size()){
            order.clear(); //in this case, rank doesnt mean anything
            /*
            for(Int32 i = 0; i<Extrema.size(); i++){
                order.push_back(i);
            }*/
            std::vector<std::string> ids = candResults->ExtremaIDs;
            Int32 ret = FixRanksUsingSortedIDs(order, ids);
            if(ret<0){
                Log.LogError("   CLineModelExtremaResult::FixRanksUsingSortedIDs failed!");
                throw std::runtime_error("   CLineModelExtremaResult::FixRanksUsingSortedIDs failed!");
            }
        }
    }
  stream << "{"<< std::endl;
  // save extrema list, on 1 line
  SaveTFloat64List(stream,"z_extrema", Extrema, order);
  stream << "," << std::endl;

  //warning: for results of second pass, no need to reorder ids cause already done in pdfcandidateresult
  //while is necessary for firstpassids!!
  // save extrema list, on 1 line
  SaveStringVector(stream,"z_extremaIDs", ExtremaIDs, order);
  stream << "," << std::endl;

  if(!zeros){
    // save extremum final rank list, on 1 line
    std::vector<int> finalRanks(Extrema.size());
    std::iota(finalRanks.begin(), finalRanks.end(), 0);
    SaveInt32Vector(stream,"z_ExtremaFinalRank", finalRanks, {});
    stream << "," << std::endl;

    // save extremaPDF list, on 1 line
    TFloat64List intgPDF = candResults->ValSumProba;
    SaveTFloat64List(stream,"z_extremaIntgPDF", intgPDF, order);
    stream << "," << std::endl;
  }
    // save extremaMerit list, on 1 line
  SaveTFloat64List(stream,"z_ExtremaMerit",ExtremaMerit, order);
  stream << "," << std::endl;

  if(!zeros){
  // save extremaMeritContinuum list, on 1 line
  SaveTFloat64List(stream,"z_ExtremaMeritContinuum",ExtremaMeritContinuum, order);
  stream << "," << std::endl;
  // save extrema Deltaz list, on 1 line
  SaveTFloat64List(stream,"z_ExtremaDeltaZ",DeltaZ, order);
  stream << "," << std::endl;
  // save extrema mTransposeM list, on 1 line
  SaveTFloat64List(stream,"z_mTransposeM",mTransposeM, order);
  stream << "," << std::endl;
  // save NDof list, on 1 line
  SaveInt32Vector(stream,"z_NDof",NDof, order);
  stream << "," << std::endl;
  // save CorrScaleMarg list, on 1 line
  SaveTFloat64List(stream,"z_CorrScaleMarg",CorrScaleMarg, order);
  stream << "," << std::endl;
  // save ExtremaLastPass list, on 1 line
  SaveTFloat64List(stream,"z_ExtremaLastPass",ExtremaLastPass, order);
  stream << "," << std::endl;
  // save lmfitPass list, on 1 line
  SaveTFloat64List(stream,"z_lmfitPass",lmfitPass, order);
  stream << "," << std::endl;
  // save snrHa list, on 1 line
  SaveTFloat64List(stream,"z_snrHa",snrHa, order);
  stream << "," << std::endl;
  // save lfHa list, on 1 line
  SaveTFloat64List(stream,"z_lfHa",lfHa, order);
  stream << "," << std::endl;

  // save snrOII list, on 1 line
  SaveTFloat64List(stream,"z_snrOII",snrOII, order);
  stream << "," << std::endl;
  // save lfOII list, on 1 line
  SaveTFloat64List(stream,"z_lfOII",lfOII, order);
  stream << "," << std::endl;
  // save bic list, on 1 line
  SaveTFloat64List(stream,"ext_BIC",bic, order);
  stream << "," << std::endl;
  // save posterior list, on 1 line
  SaveTFloat64List(stream,"ext_POSTERIOR",Posterior, order);
  stream << "," << std::endl;
  // save SigmaZ list, on 1 line
  SaveTFloat64List(stream,"ext_SigmaZ",SigmaZ, order);
  stream << "," << std::endl;
  // save LogArea list, on 1 line
  SaveTFloat64List(stream,"ext_LogArea",LogArea, order);
  stream << "," << std::endl;
  // save ContinuumIndexes list, on 1 line
  SaveTContinuumIndexListVector(stream,"ext_ContinuumIndexes",ContinuumIndexes, order);
  stream << "," << std::endl;	


  // save StrongELSNR list, on 1 line
  SaveTFloat64List(stream,"ext_StrongELSNR",StrongELSNR, order);
  stream << "," << std::endl;

  // save StrongELSNRAboveCut list, on 1 line
  SaveStringVectorOfVector(stream,"ext_StrongELSNRAboveCut",StrongELSNRAboveCut, order);
  stream << "," << std::endl;
}
  // save FittedTplName, on 1 line
  SaveStringVector(stream,"ext_FittedTplName",FittedTplName, order);
  stream << "," << std::endl;
  // save FittedTplAmplitude, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplAmplitude",FittedTplAmplitude, order);
  stream << "," << std::endl;
  // save FittedTplAmplitudeError, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplAmplitudeError",FittedTplAmplitudeError, order);
  stream << "," << std::endl;
  // save FittedTplMerit, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplMerit",FittedTplMerit, order);
  stream << "," << std::endl;
  // save FittedTplDustCoeff, on 1 line
  stream << std::setprecision(3);
  SaveTFloat64List(stream,"ext_FittedTplDustCoeff",FittedTplDustCoeff, order); 
  stream << "," << std::endl; 

  // save FittedTplMeiksinIdx, on 1 line
  SaveInt32Vector(stream,"ext_FittedTplMeiksinIdx",FittedTplMeiksinIdx, order);
  stream << "," << std::endl;
   
  // save FittedTplRedshift, on 1 line
   stream << std::setprecision(8);
  SaveTFloat64List(stream,"ext_FittedTplRedshift",FittedTplRedshift, order);  
  stream << "," << std::endl;

  // save FittedTplDtm, on 1 line
      stream << std::setprecision(8);
  SaveTFloat64List(stream,"ext_FittedTplDtm",FittedTplDtm, order);
  stream << "," << std::endl;
  // save FittedTplMtm, on 1 line
    stream << std::setprecision(8);
  SaveTFloat64List(stream,"ext_FittedTplMtm",FittedTplMtm, order);
  stream << "," << std::endl;

  if(!zeros){
  // save FittedTplLogPrior, on 1 line
    stream << std::setprecision(8);
  SaveTFloat64List(stream,"ext_FittedTplLogPrior",FittedTplLogPrior, order);
  stream << "," << std::endl;
  // save FittedTplpCoeffs, on 1 line
    stream << std::setprecision(8);
  SaveTFloat64ListOfList(stream,"ext_FittedTplpCoeffs",FittedTplpCoeffs, order);
  stream << "," << std::endl;

  // save FittedTplshapeName, on 1 line
  SaveStringVector(stream,"ext_FittedTplshapeName",FittedTplshapeName, order);
  stream << "," << std::endl;
  // save FittedTplshapeIsmCoeff, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplshapeIsmCoeff",FittedTplshapeIsmCoeff, order);
  stream << "," << std::endl;
  // save FittedTplshapeAmplitude, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplshapeAmplitude",FittedTplshapeAmplitude, order);
  stream << "," << std::endl;
  // save FittedTplshapeDtm, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplshapeDtm",FittedTplshapeDtm, order);
  stream << "," << std::endl;
  // save FittedTplshapeMtm, on 1 line
  SaveTFloat64List(stream,"ext_FittedTplshapeMtm",FittedTplshapeMtm, order);
  stream << "," << std::endl;

  // save OutsideLinesSTDFlux, on 1 line
  stream << std::scientific << std::setprecision(5);
  SaveTFloat64List(stream,"ext_OutsideLinesSTDFlux",OutsideLinesSTDFlux, order);
  stream << "," << std::endl;
  // save OutsideLinesSTDError, on 1 line
  stream << std::scientific << std::setprecision(5);
  SaveTFloat64List(stream,"ext_OutsideLinesSTDError",OutsideLinesSTDError, order);
  stream << "," << std::endl;
  }
  // save Elv, on 1 line
  stream << std::fixed << std::setprecision(1);
  SaveTFloat64List(stream,"ext_Elv",Elv, order);
  stream << "," << std::endl;
  // save Alv, on 1 line
  stream << std::fixed << std::setprecision(1);
  SaveTFloat64List(stream,"ext_Alv",Alv, order);

  stream <<  "}";

}

