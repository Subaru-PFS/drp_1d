#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <algorithm>    // std::sort
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <boost/filesystem.hpp>

#include <RedshiftLibrary/ray/lineprofile.h>

using namespace NSEpic;
using namespace std;
using namespace boost;
using namespace boost::filesystem;


CRayCatalog::CRayCatalog()
{

}

CRayCatalog::~CRayCatalog()
{

}

const CRayCatalog::TRayVector& CRayCatalog::GetList() const
{
    return m_List;
}

const CRayCatalog::TRayVector CRayCatalog::GetFilteredList(Int32 typeFilter, Int32 forceFilter) const
{
    {
        TRayVector filteredList;
        for( int i = 0; i< m_List.size(); i++ )
        {
            if( typeFilter == -1 || typeFilter == m_List[i].GetType()){
                if( forceFilter == -1 || forceFilter == m_List[i].GetForce()){
                    filteredList.push_back(m_List[i]);
                }
            }
        }
        return filteredList;
    }
}

const std::vector<CRayCatalog::TRayVector> CRayCatalog::ConvertToGroupList( TRayVector filteredList ) const
{

    std::vector<std::string> tags;
    for( int i = 0; i< filteredList.size(); i++ )
    {
        if(filteredList[i].GetGroupName() != "-1"){
            tags.push_back(filteredList[i].GetGroupName());
        }
    }

    // create the group tag set by removing duplicates
    std::sort( tags.begin(), tags.end() );
    tags.erase( std::unique( tags.begin(), tags.end() ), tags.end() );

    //get all group tags
    std::vector<TRayVector> fullList;
    for( Int32 itag = 0; itag<tags.size(); itag++){
        TRayVector taggedGroupList;
        for( int i = 0; i< filteredList.size(); i++ )
        {
            std::string group = filteredList[i].GetGroupName();
            if(group==tags[itag]){
                taggedGroupList.push_back(filteredList[i]);
            }
        }
        fullList.push_back(taggedGroupList);
    }
    //add the non grouped lines
    for( int i = 0; i< filteredList.size(); i++ )
    {
        std::string group = filteredList[i].GetGroupName();
        if(group=="-1"){
            TRayVector taggedGroupList;
            taggedGroupList.push_back(filteredList[i]);
            fullList.push_back(taggedGroupList);
        }
    }

    return fullList;
}


Bool CRayCatalog::Add( const CRay& r )
{
    TRayVector::iterator it;
    for( it = m_List.begin(); it != m_List.end(); ++it )
    {
      //TODO this should be a map with key defined as struct{name,position,type}
      // TODO this should be a map with key defined as ID, much better than struct{name,position,type} now that a ray has an ID
        // Can't add a line with a name + position + type that already exists in the list
        if( (*it).GetName() == r.GetName() && (*it).GetPosition() == r.GetPosition() && (*it).GetType() == r.GetType() )
            return false;
    }

    m_List.push_back( r );

    return true;
}

/**
 * @brief parse_asymfixed
 * Parse a ASYMFIXED_XX_YY_ZZ profile string
 */
//static void parse_asymfixed(std::string &profileString, Float64 &width, Float64 &alpha, Float64 &delta)
static void parse_asymfixed(std::string &profileString, TAsymParams &asymParams)
{
    std::vector< Float64 > numbers;
    std::string temp;
    size_t start = 0, end;

    end = profileString.find("_", start);
    if (profileString.substr(start, end-start+1) != "ASYMFIXED_") {
        throw std::runtime_error("XXX");
    }
    start = end + 1;
    end = profileString.find("_", start);
    temp = profileString.substr(start, end-start);
    asymParams.width = std::stod(temp);

    start = end + 1;
    end = profileString.find("_", start);
    temp = profileString.substr(start, end-start);
    asymParams.alpha = std::stod(temp);

    start = end + 1;
    end = profileString.length();
    temp = profileString.substr(start, end-start);
    asymParams.delta = std::stod(temp);
}

/**
 * @brief CRayCatalog::Load
 * Loads a line catalog in TSV format as of v0.4
 * @param filePath
 */
void CRayCatalog::Load( const char* filePath )
{
    std::ifstream file;

    // Clear current line list
    m_List.clear();

    if ( !exists( filePath ) ) {
      Log.LogError("Can't load line catalog : %s does not exist.", filePath);
      throw runtime_error("Can't load line catalog");
    }

    file.open( filePath, std::ifstream::in );
    if( file.rdstate() & ios_base::failbit ) {
      throw runtime_error("file cannot be opened");
    }

    string line;

    Float64 ver = -1.0;
    // Read file line by line
    while( getline( file, line ) )
    {
        // manage version
        if(ver==-1.0){
            if(line.find("#version:0.4.0")!= std::string::npos)
            {
                ver = 0.4;
            }else if(line.find("#version:0.3.0")!= std::string::npos)
            {
                ver = 0.3;
            }else if(line.find("#version:")!= std::string::npos)
            {
                Log.LogDebug("Loading Line Catalog: unable to parse version line: %s", line.c_str());
            }
            continue;
        }

        if(ver!=0.4)
        {
            Log.LogError("Line catalog version (found ver=%.3f) is not supported", ver);
            throw runtime_error("Line catalog version is not supported");
        }

        // remove comments
        if(line.compare(0,1,"#",1)==0){
            continue;
        }
        char_separator<char> sep("\t");

        // Tokenize each line
        typedef tokenizer< char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );

        // Check if it's not a comment
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() && *it != "#" )
        {
            // Parse position
            double pos = 0.0;
            try
            {
                pos = lexical_cast<double>(*it);
            }
            catch (const bad_lexical_cast& e)
            {
	      Log.LogError("Bad file format : %s [%s]", filePath, e.what());
	      throw runtime_error("Bad file format");
            }

            // Parse name
            ++it;
            string name;
            if( it != tok.end() )
            {
                name = *it;
            }
            else
            {
	      Log.LogError("Bad name in : %s", + filePath);
	      throw runtime_error("Bad name");
            }

            // Parse type
            int Etype = 0;
            ++it;
            string type = "None";
            if( it != tok.end() )
                type = *it;
            if( strcmp(type.c_str(),"A")==0 ){
                Etype = 1;
            }else if( strcmp(type.c_str(),"E")==0 ){
                Etype = 2;
            }

            // Parse weak or strong
            int Eforce = 0;
            ++it;
            string strong = "None";
            if( it != tok.end() )
                strong = *it;
            if( strcmp(strong.c_str(),"W")==0 ){
                Eforce = 1;
            }else if( strcmp(strong.c_str(),"S")==0 ){
                Eforce = 2;
            }
            else
              {
                Log.LogError("Bad force in : %s", + filePath);
                throw runtime_error("Bad force");
              }

            std::string profileName = "SYM";
            TAsymParams asymParams = {NAN, NAN, NAN};
            std::string groupName = "-1";
            Float64 nominalAmplitude = 1.0;
            std::string velGroupName = "-1";

            std::shared_ptr<CLineProfile> profile;
            // Parse profile name
            ++it;
            if( it != tok.end() ){
                profileName = *it;
                if (profileName.find("ASYMFIXED") != std::string::npos) {
                    profile = std::make_shared<CLineProfileASYMFIXED>();
                    parse_asymfixed(profileName, asymParams);
                }
                else if (profileName == "SYM") {
                    profile = std::make_shared<CLineProfileSYM>();
                }
                else if (profileName == "LOR") { 
                    profile = std::make_shared<CLineProfileLOR>();
                }
                else if (profileName == "ASYM") { 
                    profile = std::make_shared<CLineProfileASYM>();
                }
                else if (profileName == "ASYMFIT") { 
                    profile = std::make_shared<CLineProfileASYMFIT>();
                }
            }

            // Parse group name
            ++it;
            if( it != tok.end() )
            {
                groupName = *it;
                // Parse group line nominal amplitude
                ++it;
                if( it != tok.end() )
                {
                    try
                    {
                        nominalAmplitude = lexical_cast<double>(*it);
                    }
                    catch (bad_lexical_cast&)
                    {
                        Log.LogError( "Unable to read nominal amplitude value from file, setting as default (1.0)." );
                        nominalAmplitude = 1.0;
                    }
                }
                if(groupName=="" || groupName=="-1")
                {
                    nominalAmplitude = 1.0;
                }
            }

            // Parse velocity group name
            ++it;
            if( it != tok.end() ){
                velGroupName = *it;
            }

	    ++it;
	    int id=-1;
	    if( it != tok.end() ){
	      try
		{
		  id = lexical_cast<int>(*it);
		}
	      catch (const bad_lexical_cast& e)
		{
		  Log.LogError("Bad file format : %s [%s]", filePath, e.what());
		  throw runtime_error("Bad file format");
		}
            }
	    

            if( nominalAmplitude>0.0 ) //do not load a line with nominal amplitude = ZERO
            {
	      Add( CRay(name, pos, Etype, profile, Eforce, -1, -1, -1, -1, -1, -1, groupName, nominalAmplitude, velGroupName, asymParams, id) );
            }
        }
    }
    file.close();

    if(ver<0.0)
    {
        Log.LogError("Invalid line catalog file (found ver=%.3f)", ver);
        throw runtime_error("Invalid line catalog file");
    }
}

Bool CRayCatalog::Save( const char* filePath )
{
    std::ofstream file;
    file.open( filePath, std::ofstream::out );
    if( file.rdstate() & ios_base::failbit )
    {
        return false;
    }

    file << "#version:0.4.0" << std::endl;
    file << "#lambda" << "\t" << "name" << "\t" << "type" << "\t" << "force" << "\t" << "profile" << "\t" << "amp_group" << "\t" << "nominal_ampl" << "\t" << "vel_group" << std::endl;
    for( int i = 0; i< m_List.size(); i++ )
    {
        file << m_List[i].GetPosition() << "\t";
        file << m_List[i].GetName() << "\t";
        if(m_List[i].GetType() == 1)
        {
            file << "A" << "\t";
        }else
        {
            file << "E" << "\t";
        }

        if(m_List[i].GetForce() == 1)
        {
            file << "W" << "\t";
        }else
        {
            file << "S" << "\t";
        }
        file << m_List[i].GetName() << "\t";

        file << m_List[i].GetGroupName() << "\t";
        file << m_List[i].GetNominalAmplitude() << "\t";

        file << m_List[i].GetVelGroupName();

        file << std::endl;

    }
    file.close();
    return true;
}

void CRayCatalog::Sort()
{
    sort(m_List.begin(), m_List.end());
}

void CRayCatalog::ConvertVacuumToAir()
{
    TRayVector::iterator it;
    for( it = m_List.begin(); it != m_List.end(); ++it )
    {
        (*it).ConvertVacuumToAir();
    }

    return;
}
