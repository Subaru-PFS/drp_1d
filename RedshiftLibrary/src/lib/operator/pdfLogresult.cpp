#include <RedshiftLibrary/operator/pdfLogresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

using namespace NSEpic;


CPdfLogResult::CPdfLogResult()
{

}

CPdfLogResult::~CPdfLogResult()
{

}

void CPdfLogResult::SetSize( UInt32 n )
{
    Redshifts.resize(n);
    valProbaLog.resize(n);
    Overlap.resize(n);
    Status.resize(n);
}

void CPdfLogResult::Load( std::istream& stream )
{
    // Clear current list
    Redshifts.clear();
    valProbaLog.clear();
    Overlap.clear();
    Status.clear();

    std::string line;

    // Read file line by line
    while( std::getline( stream, line ) )
    {
        boost::char_separator<char> sep(" \t");

        // Tokenize each line
        typedef boost::tokenizer< boost::char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );

        // Check if it's not a comment
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() && *it != "#" )
        {
            // Parse
            Float64 r = 0.0;
            try
            {
                r = boost::lexical_cast<Float64>(*it);
            }
            catch (boost::bad_lexical_cast&)
            {
                return;
            }

            // Parse
            ++it;
            Float64 c = 0.0;
            if( it != tok.end() )
            {
                try
                {
                    c = boost::lexical_cast<Float64>(*it);
                }
                catch (boost::bad_lexical_cast&)
                {
                    return;
                }
            }
            else
            {
                return;
            }

            // Parse
            ++it;
            Float64 o = 0.0;
            if( it != tok.end() )
            {
                try
                {
                    o = boost::lexical_cast<Float64>(*it);
                }
                catch (boost::bad_lexical_cast&)
                {
                    return;
                }
            }
            else
            {
                return;
            }

            Redshifts.push_back( r );
            valProbaLog.push_back( c );
            Overlap.push_back( o );
            Status.push_back( COperator::nStatus_OK );
        }
    }
}

void CPdfLogResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream <<  "#Redshifts\tPdfLog\tOverlap"<< std::endl;
    for ( Int32 i=0; i<Redshifts.size(); i++)
    {
        stream.precision(10);
        stream
			<<  Redshifts[i] << "\t"
        		<< std::scientific << valProbaLog[i] ;
        stream.precision(5);
        stream
			 << "\t" << Overlap[i]  // << std::fixed << "\t" << Overlap[i]
			<< std::endl;
    }
}

void CPdfLogResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
	stream.precision(20);
    stream << "CPdfLogResult" << "\t" << Redshifts.size() << std::endl;
}
