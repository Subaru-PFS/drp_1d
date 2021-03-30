#include <RedshiftLibrary/operator/templatefittingresult.h>
#include <RedshiftLibrary/extremum/extremum.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <boost/algorithm/string/predicate.hpp>

using namespace NSEpic;

CTemplateFittingResult::CTemplateFittingResult()
{

}

CTemplateFittingResult::~CTemplateFittingResult()
{

}

void CTemplateFittingResult::Init(UInt32 n , Int32 nISM, Int32 nIGM)
{
    ChiSquare.resize( n );
    FitAmplitude.resize( n );
    FitAmplitudeError.resize( n );
    FitAmplitudeNegative.resize( n);
    FitDustCoeff.resize( n );
    FitMeiksinIdx.resize( n );
    FitDtM.resize( n );
    FitMtM.resize( n );
    LogPrior.resize( n );
    Redshifts.resize( n );
    Overlap.resize( n );
    Status.resize( n );

    ChiSquareIntermediate.clear();
    IsmDustCoeffIntermediate.clear();
    IgmMeiksinIdxIntermediate.clear();
    for(Int32 k=0; k<n; k++)
    {
        std::vector<TFloat64List> _ChiSquareISMList;
        std::vector<TFloat64List> _IsmDustCoeffISMList;
        std::vector<TInt32List> _IgmMeiksinIdxISMList;
        for(Int32 kism=0; kism<nISM; kism++)
        {

            TFloat64List _chi2IGMList(nIGM, DBL_MAX);
            _ChiSquareISMList.push_back(_chi2IGMList);

            TFloat64List _dustCoeffIGMList(nIGM, -1.0);
            _IsmDustCoeffISMList.push_back(_dustCoeffIGMList);

            TInt32List _meiksinIdxIGMList(nIGM, -1);
            _IgmMeiksinIdxISMList.push_back(_meiksinIdxIGMList);

        }
        ChiSquareIntermediate.push_back(_ChiSquareISMList);
        IsmDustCoeffIntermediate.push_back(_IsmDustCoeffISMList);
        IgmMeiksinIdxIntermediate.push_back(_IgmMeiksinIdxISMList);
    }
}

void CTemplateFittingResult::Load( std::istream& stream )
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

void CTemplateFittingResult::Save(std::ostream& stream ) const
{
    stream <<  "#Redshifts\tChiSquare\tOverlap"<< std::endl;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        stream <<  Redshifts[i] << std::setprecision(16) << "\t" << std::scientific << ChiSquare[i] << std::fixed << "\t" << Overlap[i] << std::endl;
    }


}

void CTemplateFittingResult::SaveLine(std::ostream& stream ) const
{
    stream << "ChisquareResult" << "\t" << Redshifts.size() << std::endl;
}
