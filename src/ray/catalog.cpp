#include <epic/redshift/ray/catalog.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>

using namespace __NS__;
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

Bool CRayCatalog::Add( const CRay& r )
{
    TRayVector::iterator it;
    for( it = m_List.begin(); it != m_List.end(); ++it )
    {
        // Can't add a ray with a name that already exists in the list
        if( (*it).GetName() == r.GetName() )
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
            ++it;
            string type = "E";
            if( it != tok.end() )
                type = *it;

            // Parse weak or strong
            ++it;
            string strong = "S";
            if( it != tok.end() )
                strong = *it;
            

            Add( CRay( name, pos, 0 ) );
        }
    }

    return true;
}
