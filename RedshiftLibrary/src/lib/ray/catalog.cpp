#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <algorithm>    // std::sort
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>

using namespace NSEpic;
using namespace std;
using namespace boost;


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
        // Can't add a line with a name + position + type that already exists in the list
        if( (*it).GetName() == r.GetName() && (*it).GetPosition() == r.GetPosition() && (*it).GetType() == r.GetType() )
            return false;
    }

    m_List.push_back( r );

    return true;
}

/**
 * @brief CRayCatalog::Load
 * Loads a line catalog in TSV format as of v0.4
 * @param filePath
 */
void CRayCatalog::Load( const char* filePath )
{
    ifstream file;

    // Clear current line list
    m_List.clear();

    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit )
      throw std::string("file cannot be opened");

    string line;

    Float64 ver = -1.0;
    // Read file line by line
    while( getline( file, line ) )
    {
        // manage version
        if(ver==-1.0){
            if(line.find("#version:0.4.0")!= std::string::npos){
                ver = 0.4;
            }
            continue;
        }

        if(ver!=0.4)
	  throw string("Line catalog version is not supported : ") + filePath;

        // remove comments
        if(line.compare(0,1,"#",1)==0){
            continue;
        }
        char_separator<char> sep(" \t");

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
            catch (bad_lexical_cast)
            {
                pos = 0.0;
                throw string("Bad file format : ") + filePath;
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
                throw string("Bad name in : ") + filePath;
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

            std::string profileName = "SYM";
            std::string groupName = "-1";
            Float64 nominalAmplitude = 1.0;
            std::string velGroupName = "-1";
            // Parse profile name
            ++it;
            if( it != tok.end() ){
                profileName = *it;
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
                    catch (bad_lexical_cast)
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

            if( nominalAmplitude>0.0 ) //do not load a line with nominal amplitude = ZERO
            {
                Add( CRay(name, pos, Etype, profileName, Eforce, -1, -1, -1, -1, -1, -1, groupName, nominalAmplitude, velGroupName) );
            }
        }
    }
    file.close();

    if(ver<0.0)
      throw string("Invalid catalog file : ") + filePath;
}

Bool CRayCatalog::Save( const char* filePath )
{
    ofstream file;
    file.open( filePath, ofstream::out );
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
        file << m_List[i].GetProfile() << "\t";

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
