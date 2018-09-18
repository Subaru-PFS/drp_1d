#include <RedshiftLibrary/operator/correlationresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

using namespace NSEpic;

CCorrelationResult::CCorrelationResult()
{

}

CCorrelationResult::~CCorrelationResult()
{

}

void CCorrelationResult::Load( std::istream& stream )
{
    // Clear current lines list
    Redshifts.clear();
    Correlation.clear();
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
            // Parse position
            Float64 r = 0.0;
            try
            {
                r = boost::lexical_cast<Float64>(*it);
            }
            catch (boost::bad_lexical_cast&)
            {
                return;
            }

            // Parse name
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

            // Parse name
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
            Correlation.push_back( c );
            Overlap.push_back( o );
            Status.push_back( COperator::nStatus_OK );
        }
    }
}

void CCorrelationResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream <<  "#Redshifts\tCorrelation\tOverlap"<< std::endl;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        stream <<  Redshifts[i] << "\t" << std::scientific << Correlation[i] << std::fixed << "\t" << Overlap[i] << std::endl;
    }

    if(Extrema.size()>0){
        stream <<  "#Extrema for z = {";
        for ( int i=0; i<Extrema.size(); i++)
        {
            stream <<  Extrema[i] << "\t";
        }
        stream << "}" << std::endl;
    }
}

void CCorrelationResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream << "tCorrelationResult" << "\t" << Redshifts.size() << std::endl;
}
