#include <epic/redshift/ray/catalog.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>

using namespace NSEpic;
using namespace std;
using namespace boost;

IMPLEMENT_MANAGED_OBJECT( CRayCatalog )

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

Bool CRayCatalog::Add( const CRay& r )
{
    TRayVector::iterator it;
    for( it = m_List.begin(); it != m_List.end(); ++it )
    {
        // Can't add a ray with a name + position + type that already exists in the list
        if( (*it).GetName() == r.GetName() && (*it).GetPosition() == r.GetPosition() && (*it).GetType() == r.GetType() )
            return false;
    }

    m_List.push_back( r );

    return true;
}

Bool CRayCatalog::Load( const char* filePath )
{
    ifstream file;

    // Clear current ray list
    m_List.clear();

    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit )
        return false;

    string line;

    // Read file line by line
    while( getline( file, line ) )
    {
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
                return false;
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
                return false;
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

            Add( CRay( name, pos, Etype, Eforce ) );
        }
    }

    return true;
}
