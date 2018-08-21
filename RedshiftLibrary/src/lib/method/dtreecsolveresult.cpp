#include <RedshiftLibrary/method/dtreecsolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>

#include <RedshiftLibrary/operator/linemodelresult.h>

#include <RedshiftLibrary/log/log.h>

#include <float.h>

using namespace NSEpic;


CDTreeCSolveResult::CDTreeCSolveResult()
{

}

CDTreeCSolveResult::~CDTreeCSolveResult()
{

}

void CDTreeCSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;
    std::string dtreepath;

    GetBestRedshift( store, redshift, merit, tplName, dtreepath );

    stream <<  "#Redshifts\tMerit\tTemplate\tMethod\tDeltaz"<< std::endl;

    stream  << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                   << "dtreec_" << dtreepath.c_str() << "\t"
                   << "-1" << "\t" //deltaz
                   << std::endl;

}

void CDTreeCSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;
    std::string dtreepath;

    GetBestRedshift( store, redshift, merit, tplName, dtreepath );

    stream  << store.GetSpectrumName() << "\t"
            << store.GetProcessingID() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "dtreec_" << dtreepath.c_str() << "\t"
                << "-1" << "\t" //deltaz
                << "-1" << std::endl; //reliability label



    Log.LogInfo( "DecisionalTreeC Solution: best z found = %.5f", redshift);
}


Bool CDTreeCSolveResult::GetBestRedshift(const CDataStore& store, Float64& redshift, Float64& merit , std::string &tplName, std::string &dtreepath) const
{
    //*
    // combined merit curve minimum
    std::string scope_lincomb = store.GetScope( *this ) + "dtreeCsolve.resultdtreeCCombined";

    auto _results = store.GetGlobalResult(scope_lincomb.c_str());
    auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>(_results.lock());



    //estimate the combined chisquare result
    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = -1.0;
    for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
    {
        Float64 post = meritResult->ChiSquare[i];
        if( post < tmpMerit )
        {
            tmpMerit = post;
            tmpRedshift = meritResult->Redshifts[i];
            tplName = GetBestContinuumTplNameAtRedshift( store, tmpRedshift);
        }
    }


    redshift = tmpRedshift;
    merit = tmpMerit;
    dtreepath = "0.0";
    return true;
    //*/
}


std::string CDTreeCSolveResult::GetBestContinuumTplNameAtRedshift( const CDataStore& store, Float64 z) const
{
    // combined merit curve minimum
    std::string scope = store.GetScope( *this ) + "dtreeCsolve.chisquare2solve." + m_chi2ScopeStr.c_str();


    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    std::string tplName = "";


    //find best merit for this z
    Float64 iz = -1;
    Float64 tmpMerit = DBL_MAX ;
    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>((*it).second);

        iz = -1;
        if(iz<0)
        {
            Float64 minzdiff = DBL_MAX;
            for( Int32 i=0; i<meritResult->Redshifts.size(); i++ )
            {
                Float64 diff = fabs(meritResult->Redshifts[i]-z);
                if(minzdiff>diff)
                {
                    minzdiff = diff;
                    iz = i;
                }
            }
        }

        if( tmpMerit > meritResult->ChiSquare[iz]){
            tmpMerit = meritResult->ChiSquare[iz];
            tplName = (*it).first;
        }

    }

    return tplName;

}
