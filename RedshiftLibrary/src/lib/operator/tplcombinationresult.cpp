#include <RedshiftLibrary/operator/tplcombinationresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <boost/algorithm/string/predicate.hpp>

using namespace NSEpic;

CTplcombinationResult::CTplcombinationResult()
{

}

CTplcombinationResult::~CTplcombinationResult()
{

}

void CTplcombinationResult::Init(UInt32 n , Int32 nISM, Int32 nIGM)
{
    ChiSquare.resize( n );
    FitAmplitude.resize( n );
    FitDustCoeff.resize( n );
    FitMeiksinIdx.resize( n );
    FitDtM.resize( n );
    FitMtM.resize( n );
    Redshifts.resize( n );
    Overlap.resize( n );
    Status.resize( n );

    for(Int32 k=0; k<n; k++)
    {
        std::vector<TFloat64List> _ChiSquareISMList;
        for(Int32 kism=0; kism<nISM; kism++)
        {

            TFloat64List _chi2List(nIGM, DBL_MAX);
            _ChiSquareISMList.push_back(_chi2List);

        }
        ChiSquareIntermediate.push_back(_ChiSquareISMList);
    }
}

void CTplcombinationResult::Load( std::istream& stream )
{
    // Clear current lines list
    Redshifts.clear();
    ChiSquare.clear();
    Overlap.clear();
    Status.clear();

    std::string line;

    // Read file line by line
    while( std::getline( stream, line ) )
    {
        if( !boost::starts_with( line, "#" ) )
        {
            boost::char_separator<char> sep(" \t");

            // Tokenize each line
            typedef boost::tokenizer< boost::char_separator<char> > ttokenizer;
            ttokenizer tok( line, sep );

            // Check if it's not a comment
            ttokenizer::iterator it = tok.begin();
            if( it != tok.end() )
            {
                // Parse redshift
                Float64 r = -1.0;
                try
                {
                    r = boost::lexical_cast<Float64>(*it);
                }
                catch (boost::bad_lexical_cast&)
                {
                    return;
                }

                // Parse merit
                ++it;
                Float64 c = DBL_MAX;
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

                // Parse overlap
                ++it;
                Float64 o = -1.0;
                if( it != tok.end() )
                {
                    try
                    {
                        o = boost::lexical_cast<Float64>(*it);
                    }
                    catch (boost::bad_lexical_cast&)
                    {
                        o = -1.0;
                    }
                }
                else
                {
                    o = -1.0;
                }

                Redshifts.push_back( r );
                ChiSquare.push_back( c );
                Overlap.push_back( o );
                Status.push_back( COperator::nStatus_OK );
            }
        }
    }
}

void CTplcombinationResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream <<  "#Redshifts\tChiSquare\tOverlap"<< std::endl;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        stream <<  Redshifts[i] << std::setprecision(16) << "\t" << std::scientific << ChiSquare[i] << std::fixed << "\t" << Overlap[i] << std::endl;
    }

    if(Extrema.size()>0){
        stream <<  "#Extrema for z = {";
        for ( int i=0; i<Extrema.size(); i++)
        {
            stream <<  Extrema[i] << "\t";
        }
        stream << "}" << std::endl;

        if(ChiSquare.size()>0){
            stream <<  "#ExtremaMerit = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << std::setprecision(16) << "\t" << std::scientific <<   ChiSquare[idx] << "\t";
            }
            stream << "}" << std::endl;
        }


        if(FitAmplitude.size()>0){
            stream <<  "#Extrema FitAmplitudes = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << std::setprecision(16) << "\t" << std::scientific <<   FitAmplitude[idx] << "\t";
            }
            stream << "}" << std::endl;
        }

        if(FitDustCoeff.size()>0){
            stream <<  "#Extrema FitDustCoeff = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << std::setprecision(4) << "\t" <<   FitDustCoeff[idx] << "\t";
            }
            stream << "}" << std::endl;
        }

        if(FitMeiksinIdx.size()>0){
            stream <<  "#Extrema FitMeiksinIdx = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << "\t" <<   FitMeiksinIdx[idx] << "\t";
            }
            stream << "}" << std::endl;
        }
    }
}

void CTplcombinationResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream << "ChisquareResult" << "\t" << Redshifts.size() << std::endl;
}
