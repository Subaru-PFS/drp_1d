#include <RedshiftLibrary/processflow/resultstore.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/template.h>

#include <boost/filesystem.hpp>
#include <fstream>
#include <boost/algorithm/string.hpp>


using namespace NSEpic;

namespace bfs = boost::filesystem;


COperatorResultStore::COperatorResultStore()
{

}

COperatorResultStore::~COperatorResultStore()
{

}

void COperatorResultStore::StoreResult( TResultsMap& map, const std::string& path, const std::string& name,
                                        std::shared_ptr<const COperatorResult> result )
{
    std::string scopedName;
    if( ! path.empty() ) {
        scopedName = path;
        scopedName.append( "." );
    }
    scopedName.append( name );

    TResultsMap::iterator it = map.find( name );
    if( it != map.end() )
    {
        DebugError( "Result already exist");
        return;
    }

    map[ scopedName ] = result;
}

void COperatorResultStore::StorePerTemplateResult( const CTemplate& t, const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result )
{
    TPerTemplateResultsMap::iterator it = m_PerTemplateResults.find( t.GetName() );
    if( it == m_PerTemplateResults.end() )
    {
        m_PerTemplateResults[ t.GetName() ] = TResultsMap();
    }

    StoreResult( m_PerTemplateResults[ t.GetName() ], path, name, result );
}

void COperatorResultStore::StoreGlobalResult( const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result )
{
    StoreResult( m_GlobalResults, path, name, result );
}

std::weak_ptr<const COperatorResult> COperatorResultStore::GetPerTemplateResult( const CTemplate& t, const std::string& name ) const
{
    TPerTemplateResultsMap::const_iterator it1 = m_PerTemplateResults.find( t.GetName() );
    if( it1 != m_PerTemplateResults.end() )
    {
        const TResultsMap& m = (*it1).second;
        TResultsMap::const_iterator it2 = m.find( name );
        if( it2 != m.end() )
        {
            return (*it2).second;
        }
    }

    return std::weak_ptr<const COperatorResult>();
}

TOperatorResultMap COperatorResultStore::GetPerTemplateResult( const std::string& name ) const
{
    TOperatorResultMap map;
    TPerTemplateResultsMap::const_iterator it;
    for( it = m_PerTemplateResults.begin(); it != m_PerTemplateResults.end(); ++it )
    {
        std::string tplName = (*it).first;

        const TResultsMap& m = (*it).second;
        TResultsMap::const_iterator it2 = m.find( name );
        if( it2 != m.end() )
        {
            map[tplName] = (*it2).second;
        }
    }

    return map;
}

/**
 * /brief Returns the best global result, if there is one.
 * 
 * The somewhat strange syntax of this method is due to the usage of std::map, which will yield an iterator when accessing a member using .find().
 * This method will find the global result entry with key "name", and if it exists, will return its .second member (which will be a COperatorResult or derived object). Otherwise, will return a pointer to an empty COperatorResult.
 */

std::weak_ptr<const COperatorResult> COperatorResultStore::GetGlobalResult( const std::string& name ) const
{
    TResultsMap::const_iterator it = m_GlobalResults.find( name );
    if( it != m_GlobalResults.end() )
    {
        return (*it).second;
    }

    return std::weak_ptr<const COperatorResult>();
}

/**
 * @brief COperatorResultStore::CreateResultStorage
 *
 * @return 0 if already exists, 1 if just created, -1 if error
 */
Int32 COperatorResultStore::CreateResultStorage( std::fstream& stream, const bfs::path& path, const bfs::path& baseDir ) const
{
    Int32 ret = -1;

    bfs::path outputFilePath = bfs::path( baseDir );
    outputFilePath /= path.string();


    if( bfs::exists( outputFilePath.parent_path() ) == false )
    {
        bfs::create_directories( outputFilePath.parent_path() );
    }

    if( bfs::exists( outputFilePath ) == false )
    {
        ret = 1;
    }else{
        ret = 0;
    }

    stream.open( outputFilePath.string().c_str(), std::fstream::out | std::fstream::app);
    if( stream.rdstate() & std::ios_base::failbit )
    {
        return -1;
    }

    return ret;
}

void COperatorResultStore::SaveRedshiftResult( const CDataStore& store, const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "redshift.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift\tMerit\tTemplate\tMethod\tDeltaz\tReliability\tsnrHa\tlfHa\tType"<< std::endl;
        }
        //*/

        auto  result = GetGlobalResult( "redshiftresult" ).lock();
        if(result){
            result->SaveLine( store, outputStream );
        }else{
            throw std::runtime_error("Unable to retrieve redshift result for saving");
        }
    }
}

void COperatorResultStore::SaveRedshiftResultError(  const std::string spcName,  const std::string processingID, const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "redshift.csv" ), dir );

        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift\tMerit\tTemplate\tMethod\tDeltaz\tReliability\tsnrHa\tlfHa\tType"<< std::endl;
        }


        outputStream <<  spcName << "\t" << processingID << "\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1"<< std::endl;
    }
}

void COperatorResultStore::SaveStellarResult( const CDataStore& store, const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "stellar.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift\tMerit\tTemplate\tMethod\tDeltaz\tReliability\tsnrHa\tlfHa\tType"<< std::endl;
        }
        //*/

        auto  result = GetGlobalResult( "stellarsolve.stellarresult" ).lock();
        if(result){
            result->SaveLine( store, outputStream );
        }//else{
        //    throw std::runtime_error("Unable to retrieve stellar result for saving");
        //}
    }
}

void COperatorResultStore::SaveStellarResultError(  const std::string spcName,  const std::string processingID, const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "stellar.csv" ), dir );

        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift\tMerit\tTemplate\tMethod\tDeltaz\tReliability\tsnrHa\tlfHa\tType"<< std::endl;
        }


        outputStream <<  spcName << "\t" << processingID << "\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1"<< std::endl;
    }
}


void COperatorResultStore::SaveClassificationResult( const CDataStore& store, const bfs::path& dir )
{
    // Append classif. result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "classification.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\ttype\tevidenceG\tevidenceS\tevidenceQ"<< std::endl;
        }
        //*/

        auto  result = GetGlobalResult( "classificationresult" ).lock();
        if(result){
            result->SaveLine( store, outputStream );
        }//else{
        //    throw std::runtime_error("Unable to retrieve classification result for saving");
        //}
    }
}

void COperatorResultStore::SaveClassificationResultError(  const std::string spcName,  const std::string processingID, const bfs::path& dir )
{
    // Append best classification result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "classification.csv" ), dir );

        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\ttype\tevidenceG\tevidenceS\tevidenceQ"<< std::endl;
        }


        outputStream <<  spcName << "\t" << processingID << "\t-1\t-1\t-1\t-1"<< std::endl;
    }
}

void COperatorResultStore::SaveCandidatesResult( const CDataStore& store, const bfs::path& dir )
{
    // Append candidate result line to output candidate file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "candidates.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRank_1\tRedshift_1\tProb_1\tgaussAmp_1\tgaussSigma_1\tRank_2\tRedshift_2\tProb_2\tgaussAmp_2\tgaussSigma_2\t..."<< std::endl;
        }
        //*/

        auto  result = GetGlobalResult( "candidatesresult" ).lock();
        if(result){
            result->SaveLine( store, outputStream );
        }
    }
}

void COperatorResultStore::SaveCandidatesResultError( const std::string spcName,  const std::string processingID, const bfs::path& dir )
{
    // Append candidate result line to output candidate file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "candidates.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift_1\tProb_1\tgaussAmp_1\tgaussSigma_1\tRedshift_2\tProb_2\tgaussAmp_2\tgaussSigma_2\t..."<< std::endl;
        }
        //*/

        outputStream <<  spcName << "\t" << processingID << "\t-1\t-1\t-1\t-1\t-1\t-1\t-1"<< std::endl;
    }
}

void COperatorResultStore::SaveReliabilityResult( const CDataStore& store, const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        CreateResultStorage( outputStream, bfs::path( "redshiftreliability.csv" ), dir );

        auto  result = GetGlobalResult( "zReliability/result.zpredict" ).lock();
        if(result){
            result->SaveLine( store, outputStream );
        }
    }
}


void COperatorResultStore::SaveAllResults( const CDataStore& store, const bfs::path& dir, const std::string opt ) const
{
    std::string opt_lower = opt;
    boost::algorithm::to_lower(opt_lower);

    // Store global result
    {
        TResultsMap::const_iterator it;
        for( it=m_GlobalResults.begin(); it != m_GlobalResults.end(); it++ )
        {
            std::string resultName = (*it).first;
            auto  result = (*it).second;

            bool saveThisResult = false;
            if(opt_lower=="all" || opt_lower=="global"){
                saveThisResult = true;
            }else if(opt_lower=="linemeas")
            {
                std::string linemeasTagRes = "linemodel_fit_extrema_0";
                std::size_t foundstr = resultName.find(linemeasTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }
            }else if(opt_lower=="default")
            {
                //save best fit parameters
                std::string linemeasTagRes = "linemodel_fit_extrema_0";
                std::size_t foundstr = resultName.find(linemeasTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }
                //save best model
                std::string linemodelTagRes = "linemodel_spc_extrema_0";
                foundstr = resultName.find(linemodelTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }
                //save best continuum
                std::string continuummodelTagRes = "linemodel_continuum_extrema_0";
                foundstr = resultName.find(continuummodelTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }
                //save redshiftresult brief
                std::string redshiftresultTagRes = "redshiftresult";
                foundstr = resultName.find(redshiftresultTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }

                //save pdf
                std::string pdfTagRes = "logposterior.logMargP_Z_data";
                foundstr = resultName.find(pdfTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }
            }

            if(!saveThisResult)
            {
                continue;
            }

            std::fstream outputStream;
            // Save result at root of output directory
            CreateResultStorage( outputStream, bfs::path( resultName + ".csv"), dir );
            result->Save( store, outputStream );
        }
    }

    // Store per template results
    if(opt_lower=="all"){
        TPerTemplateResultsMap::const_iterator it;
        for( it=m_PerTemplateResults.begin(); it != m_PerTemplateResults.end(); it++ )
        {
            std::string templateName = (*it).first;
            const TResultsMap& resultMap = (*it).second;

            TResultsMap::const_iterator it2;
            for( it2=resultMap.begin(); it2 != resultMap.end(); it2++ )
            {
                std::string resultName = (*it2).first;
                auto  result = (*it2).second;

                std::fstream outputStream;
                // Save result in sub directories of output directory
                bfs::path outputFilePath = bfs::path( templateName );
                outputFilePath /= std::string( resultName + ".csv" );
                CreateResultStorage( outputStream, outputFilePath, dir );
                result->Save( store, outputStream );
            }
        }
    }
}


std::string COperatorResultStore::GetScope( const COperatorResult&  result) const
{
    std::string n="";

    TResultsMap::const_iterator it;
    for( it=m_GlobalResults.begin(); it != m_GlobalResults.end(); it++ )
    {
        auto  r = (*it).second;
        if(&result == r.get() ){
            std::string s = (*it).first;

            std::size_t found = s.rfind(".");
            if (found!=std::string::npos){
                n = s.substr(0,found);
                n = n.append(".");
            }
            break;
        }
    }

    return n;
}
