#include <RedshiftLibrary/method/chisquarelogsolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <stdio.h>
#include <float.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>


using namespace NSEpic;

CChisquareLogSolveResult::CChisquareLogSolveResult()
{
    m_type = nType_raw;
}

CChisquareLogSolveResult::~CChisquareLogSolveResult()
{

}

Void CChisquareLogSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;
    Float64 amp;
    Float64 dustCoeff;
    Int32 meiksinIdx;

    GetBestRedshift( store, redshift, merit, tplName, amp, dustCoeff, meiksinIdx );

    stream <<  "#Redshifts\tMerit\tTemplate\tampl\tdustcoeff\tmeiksinidx"<< std::endl;

    stream  << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << amp << "\t"
                << dustCoeff << "\t"
                << meiksinIdx << std::endl;


    stream << std::endl;
    stream << std::endl;
    std::string detailStr;
    GetBestRedshiftPerTemplateString( store, detailStr);

    stream << detailStr.c_str();
}

Bool CChisquareLogSolveResult::GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const
{
    std::string scopeStr;
    if(m_type == nType_raw){
        scopeStr = "chisquare";
    }else if(m_type == nType_all){
        scopeStr = "chisquare";
    }else if(m_type == nType_noContinuum){
        scopeStr = "chisquare_nocontinuum";
    }else if(m_type == nType_continuumOnly){
        scopeStr = "chisquare_continuum";
    }
    std::string scope = store.GetScope( *this ) + "chisquarelogsolve." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());


    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        Float64 tmpMerit = DBL_MAX ;
        Float64 tmpRedshift = 0.0;
        std::string tmpTplName;

        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>( (*it).second );
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritResult->ChiSquare[i] < tmpMerit && meritResult->Status[i] == COperator::nStatus_OK )
            {
                tmpMerit = meritResult->ChiSquare[i];
                tmpRedshift = meritResult->Redshifts[i];
                tmpTplName = (*it).first;
            }
        }



        if( tmpMerit < DBL_MAX )
        {
            char tmpChar[256];
            sprintf(tmpChar, "%f\t%f\t%s\n", tmpRedshift, tmpMerit, tmpTplName.c_str());
            output.append(tmpChar);
        }else{
            char tmpChar[256];
            sprintf(tmpChar, "-1\t-1\t%s\n", tmpTplName.c_str());
            output.append(tmpChar);
        }
    }


    return true;

}

Void CChisquareLogSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    char tmpChar[256];
    Float64 dtreepathnum;
    store.GetParam( "dtreepathnum", dtreepathnum );
    sprintf(tmpChar, "%.2f", dtreepathnum);

    Float64 redshift;
    Float64 merit;
    std::string tplName="-1";

    //unused
    Float64 amplitude;
    Float64 dustCoeff;
    Int32 meiksinIdx;

    if(m_bestRedshiftMethod==0)
    {
        GetBestRedshift( store, redshift, merit, tplName, amplitude, dustCoeff, meiksinIdx );
        Log.LogInfo( "Chisquarelogsolve-result: extracting best redshift from chi2 extrema: z=%f", redshift);
    }else if(m_bestRedshiftMethod==2)
    {
        GetBestRedshiftFromPdf( store, redshift, merit );
        Log.LogInfo( "Chisquarelogsolve-result: extracting best redshift from PDF: z=%f", redshift);
        GetBestModel(store, redshift, tplName);
        Log.LogInfo( "Chisquarelogsolve-result: extracted best model: model=%s", tplName.c_str());
    }else{
        Log.LogError( "Chisquarelogsolve-result: can't parse best redshift estimation method");
    }

    stream  << store.GetSpectrumName() << "\t"
            << store.GetProcessingID() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "ChisquareLogSolve_" << tmpChar << "\t"
                << "-1" << "\t" //deltaz
                << m_ReliabilityLabel << std::endl;

}

Bool CChisquareLogSolveResult::GetBestRedshift(const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName, Float64 &amplitude, Float64 &dustCoeff, Int32 &meiksinIdx ) const
{
    std::string scopeStr;
    if(m_type == nType_raw){
        scopeStr = "chisquare";
    }else if(m_type == nType_all){
        scopeStr = "chisquare";
    }else if(m_type == nType_noContinuum){
        scopeStr = "chisquare_nocontinuum";
    }else if(m_type == nType_continuumOnly){
        scopeStr = "chisquare_continuum";
    }
    std::string scope = store.GetScope( *this ) + "chisquarelogsolve." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;
    std::string tmpTplName = "-1";
    Float64 tmpAmp = 0.0;
    Float64 tmpDustCoeff = 0.0;
    Int32 tmpMeiksinIdx = 0;

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>( (*it).second );
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritResult->ChiSquare[i] < tmpMerit && meritResult->Status[i] == COperator::nStatus_OK )
            {
                tmpMerit = meritResult->ChiSquare[i];
                tmpRedshift = meritResult->Redshifts[i];
                tmpAmp = meritResult->FitAmplitude[i];
                tmpDustCoeff = meritResult->FitDustCoeff[i];
                tmpMeiksinIdx = meritResult->FitMeiksinIdx[i];
                tmpTplName = (*it).first;
            }
        }
    }


    if( tmpMerit < DBL_MAX )
    {
        redshift = tmpRedshift;
        merit = tmpMerit;
        tplName = tmpTplName;
        amplitude = tmpAmp;
        dustCoeff = tmpDustCoeff;
        meiksinIdx = tmpMeiksinIdx;
        return true;
    }

    return false;

}


/**
 * \brief Searches the best_z = argmax(pdf)
 * output: redshift = argmax(pdf)
 * output: merit = chi2(redshift)
 *
 **/
Bool CChisquareLogSolveResult::GetBestRedshiftFromPdf( const CDataStore& store, Float64& redshift, Float64& merit ) const
{
    std::string scope_res = "zPDF/logposterior.logMargP_Z_data";
    auto results_pdf =  store.GetGlobalResult( scope_res.c_str() );
    std::shared_ptr<const CPdfMargZLogResult> logzpdf1d = std::dynamic_pointer_cast<const CPdfMargZLogResult>( results_pdf.lock() );

    if(!logzpdf1d)
    {
        Log.LogError( "GetBestRedshiftFromPdf: no pdf results retrieved from scope: %s", scope_res.c_str());
        return false;
    }


    Float64 tmpProbaLog = -DBL_MAX;
    Float64 tmpRedshift = 0.0;

    for( Int32 i=0; i<logzpdf1d->Redshifts.size(); i++ )
    {

        Float64 probaLog = logzpdf1d->valProbaLog[i];
        if(probaLog>tmpProbaLog)
        {
            tmpProbaLog = probaLog;
            tmpRedshift = logzpdf1d->Redshifts[i];
        }
    }


    redshift = tmpRedshift;
    merit = tmpProbaLog;
    return true;
}

/**
 * @brief CChisquareLogSolveResult::GetBestModel
 * Find the best tpl-name corresponding to the input arg. z and to the minimum Chisquare
 * @param store
 * @param z
 * @param tplName
 * @return
 */
Int32 CChisquareLogSolveResult::GetBestModel(const CDataStore& store, Float64 z, std::string& tplName) const
{
    tplName = "-1";
    Bool foundRedshiftAtLeastOnce = false;

    std::string scopeStr;
    if(m_type == nType_raw){
        scopeStr = "chisquare";
    }else if(m_type == nType_all){
        scopeStr = "chisquare";
    }else if(m_type == nType_noContinuum){
        scopeStr = "chisquare_nocontinuum";
    }else if(m_type == nType_continuumOnly){
        scopeStr = "chisquare_continuum";
    }
    std::string scope = store.GetScope( *this ) + "chisquarelogsolve." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    Float64 tmpMerit = DBL_MAX ;
    std::string tmpTplName = "-1";
    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>( (*it).second );
        if(meritResult->ChiSquare.size()<1){
            continue;
        }

        Int32 idx=-1;
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritResult->Status[i] == COperator::nStatus_OK )
            {
                Float64 tmpRedshift = meritResult->Redshifts[i];
                if(tmpRedshift == z){
                    idx = i;
                    foundRedshiftAtLeastOnce = true;
                    break;
                }
            }
        }
        if(idx>-1 && meritResult->ChiSquare[idx]<tmpMerit){
            tmpMerit = meritResult->ChiSquare[idx];
            tmpTplName = (*it).first;
        }

    }

    if(foundRedshiftAtLeastOnce){
        tplName = tmpTplName;
    }else{
        return -1;
    }
    return 1;
}
