#include <epic/redshift/operator/linemodelresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CLineModelResult )

CLineModelResult::CLineModelResult()
{

}

CLineModelResult::~CLineModelResult()
{

}


Void CLineModelResult::Load( std::istream& stream )
{
    // Clear current ray list
    Redshifts.clear();
    ChiSquare.clear();
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
            catch (boost::bad_lexical_cast)
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
                catch (boost::bad_lexical_cast)
                {
                    return;
                }
            }
            else
            {
                return;
            }


            Redshifts.push_back( r );
            ChiSquare.push_back( c );
            Status.push_back( COperator::nStatus_OK );
        }
    }
}

Void CLineModelResult::Save( const COperatorResultStore& store, std::ostream& stream ) const
{
    stream <<  "#Redshifts\tChiSquare\tOverlap"<< std::endl;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        stream <<  Redshifts[i] << std::setprecision(16) << "\t" << std::scientific << ChiSquare[i] << std::fixed << std::endl;
    }

    // save extrema list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#Extrema for z = {";
        for ( int i=0; i<Extrema.size(); i++)
        {
            stream <<  Extrema[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save SigmaZ list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#SigmaZ for each extrema = {";
        for ( int i=0; i<SigmaZ.size(); i++)
        {
            stream <<  SigmaZ[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save LogArea list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#LogArea for each extrema = {";
        for ( int i=0; i<LogArea.size(); i++)
        {
            stream <<  LogArea[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save linemodel solution
    if(LineModelSolutions.size()>0){
        for ( UInt32 i=0; i<Extrema.size(); i++)
        {
            stream <<  "#linemodel solution " << i << " for z = " <<  std::fixed <<  Extrema[i];
            if(LogArea.size()>i){
                stream <<  ", LogArea = " <<  LogArea[i];
            }
            Int32 idx=0;
            for ( UInt32 i2=0; i2<LineModelSolutions.size(); i2++)
            {
                if(Redshifts[i2] == Extrema[i]){
                    idx = i2;
                    break;
                }
            }
            stream << ", merit = " <<  ChiSquare[idx] << "{" <<  std::endl;
            for ( UInt32 j=0; j<LineModelSolutions[idx].Amplitudes.size(); j++)
            {
                stream <<  "#";
                std::string typeStr="";
                if(restRayList[j].GetType() == CRay::nType_Absorption){
                    typeStr = "A";
                }else{
                    typeStr = "E";
                }
                stream <<  typeStr << "\t";
                stream <<  std::fixed << restRayList[j].GetPosition() << "\t";
                stream << std::scientific <<  LineModelSolutions[idx].Amplitudes[j] << std::endl;
            }
            stream << "#}" << std::endl;
        }
    }
}

Void CLineModelResult::SaveLine( const COperatorResultStore& store, std::ostream& stream ) const
{
    stream << "LineModelResult" << "\t" << Redshifts.size() << std::endl;
}

