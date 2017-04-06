#include <epic/redshift/method/chisquare2solveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/operator/correlationresult.h>
#include <stdio.h>
#include <float.h>

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

    GetBestRedshift( store, redshift, merit, tplName );

    stream <<  "#Redshifts\tMerit\tTemplate"<< std::endl;

    stream  << redshift << "\t"
                << merit << "\t"
                << tplName << std::endl;


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
    std::string tplName;

    GetBestRedshift( store, redshift, merit, tplName );

    stream  << store.GetSpectrumName() << "\t"
            << store.GetProcessingID() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "Chisquare2Solve_" << tmpChar << "\t"
                << "-1" << "\t" //deltaz
                << std::endl;

}

Bool CChisquare2SolveResult::GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const
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
    std::string tmpTplName;

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
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
    }


    if( tmpMerit < DBL_MAX )
    {
        redshift = tmpRedshift;
        merit = tmpMerit;
        tplName = tmpTplName;
        return true;
    }

    return false;

}



