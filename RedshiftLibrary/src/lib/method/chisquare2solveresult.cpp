#include <RedshiftLibrary/method/chisquare2solveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/correlationresult.h>
#include <stdio.h>
#include <float.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>

using namespace NSEpic;

CChisquare2SolveResult::CChisquare2SolveResult()
{
    m_type = nType_raw;
}

CChisquare2SolveResult::~CChisquare2SolveResult()
{

}

Void CChisquare2SolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;
    Float64 dustCoeff;
    Int32 meiksinIdx;

    GetBestRedshift( store, redshift, merit, tplName, dustCoeff, meiksinIdx );

    stream <<  "#Redshifts\tMerit\tTemplate\tdustcoeff\tmeiksinidx"<< std::endl;

    stream  << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << std::setprecision(4) << dustCoeff << "\t"
                << meiksinIdx << std::endl;

    stream << std::endl;
    stream << std::endl;
    std::string detailStr;
    GetBestRedshiftPerTemplateString( store, detailStr);

    stream << detailStr.c_str();
}

Bool CChisquare2SolveResult::GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const
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
    std::string scope = store.GetScope( *this ) + "chisquare2solve." + scopeStr.c_str();
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

Void CChisquare2SolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    char tmpChar[256];
    Float64 dtreepathnum;
    store.GetParam( "dtreepathnum", dtreepathnum );
    sprintf(tmpChar, "%.2f", dtreepathnum);

    Float64 redshift;
    Float64 merit;
    std::string tplName="-1";

    //unused
    Float64 dustCoeff;
    Int32 meiksinIdx;

    if(m_bestRedshiftMethod==0)
    {
        GetBestRedshift( store, redshift, merit, tplName, dustCoeff, meiksinIdx );
        Log.LogInfo( "Chisquare2solve-result: extracting best redshift from chi2 extrema: z=%f", redshift);
    }else if(m_bestRedshiftMethod==2)
    {
        GetBestRedshiftFromPdf( store, redshift, merit );
        Log.LogInfo( "Chisquare2solve-result: extracting best redshift from PDF: z=%f", redshift);
    }else{
        Log.LogError( "Chisquare2solve-result: can't parse best redshift estimation method");
    }


    stream  << store.GetSpectrumName() << "\t"
            << store.GetProcessingID() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "Chisquare2Solve_" << tmpChar << "\t"
                << "-1" << "\t" //deltaz
                << m_ReliabilityLabel << std::endl; //reliability label

}

Bool CChisquare2SolveResult::GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName, Float64& dustCoeff, Int32& meiksinIdx ) const
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
    std::string scope = store.GetScope( *this ) + "chisquare2solve." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());


    Int32 maxIndex = 0;
    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;
    std::string tmpTplName = "-1";
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
Bool CChisquare2SolveResult::GetBestRedshiftFromPdf( const CDataStore& store, Float64& redshift, Float64& merit ) const
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

