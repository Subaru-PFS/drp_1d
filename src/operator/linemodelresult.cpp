#include <epic/redshift/operator/linemodelresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

using namespace NSEpic;


CLineModelResult::CLineModelResult()
{

}

CLineModelResult::~CLineModelResult()
{

}


Void CLineModelResult::Load( std::istream& stream )
{
    // Clear current lines list
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

Void CLineModelResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream <<  "#Redshifts\tChiSquare\tOverlap"<< std::endl;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        stream <<  Redshifts[i] << std::setprecision(32) << "\t" << std::scientific << ChiSquare[i] << std::fixed << std::endl;
    }

    // save extrema list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#Extrema for z = {";
        for ( int i=0; i<Extrema.size(); i++)
        {
//            if(!IsLocalExtrema[i]){
//                continue;
//            }
            stream <<  Extrema[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save bic list, on 1 line
    if(bic.size()>0){
        stream <<  "#BIC for each extrema = {";
        for ( int i=0; i<bic.size(); i++)
        {
//            if(!IsLocalExtrema[i]){
//                continue;
//            }
            stream <<  bic[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save posterior list, on 1 line
    if(Posterior.size()>0){
        stream <<  "#POSTERIOR for each extrema = {";
        for ( int i=0; i<Posterior.size(); i++)
        {
//            if(!IsLocalExtrema[i]){
//                continue;
//            }
            stream <<  Posterior[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save SigmaZ list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#SigmaZ for each extrema = {";
        for ( int i=0; i<SigmaZ.size(); i++)
        {
//            if(!IsLocalExtrema[i]){
//                continue;
//            }
            stream <<  SigmaZ[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save LogArea list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#LogArea for each extrema = {";
        for ( int i=0; i<LogArea.size(); i++)
        {
//            if(!IsLocalExtrema[i]){
//                continue;
//            }
            stream <<  LogArea[i] << "\t";
        }
        stream << "}" << std::endl;
    }
}

Void CLineModelResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream << "LineModelResult" << "\t" << Redshifts.size() << std::endl;
}

Int32 CLineModelResult::GetNLinesOverCutThreshold(Int32 extremaIdx, Float64 snrThres, Float64 fitThres) const
{
    Int32 nSol=0;
    if(Extrema.size()>extremaIdx)
    {
        Int32 solutionIdx=0;
        for ( UInt32 i2=0; i2<LineModelSolutions.size(); i2++)
        {
            if(Redshifts[i2] == Extrema[extremaIdx]){
                solutionIdx = i2;
                break;
            }
        }

        std::vector<Int32> indexesSols;
        for ( UInt32 j=0; j<LineModelSolutions[solutionIdx].Amplitudes.size(); j++)
        {
            //skip if already sol
            bool alreadysol = false;
            for(Int32 i=0; i<indexesSols.size(); i++){
                if(LineModelSolutions[solutionIdx].ElementId[j] == indexesSols[i]){
                    alreadysol=true;
                    continue;
                }
            }
            if(alreadysol){
                continue;
            }
            if(!LineModelSolutions[solutionIdx].Rays[j].GetIsStrong()){
                continue;
            }

            Float64 noise = LineModelSolutions[solutionIdx].Errors[j];
            if(noise>0){
                Float64 snr = LineModelSolutions[solutionIdx].Amplitudes[j]/noise;
                Float64 Fittingsnr = LineModelSolutions[solutionIdx].Amplitudes[j]/LineModelSolutions[solutionIdx].FittingError[j];
                if(snr>=snrThres && Fittingsnr>=fitThres){
                    nSol++;
                    indexesSols.push_back(LineModelSolutions[solutionIdx].ElementId[j]);
                }
            }
            /*
            else{
            //WARNING: this is a quick fix to deal with the case when errors are set to 0 by the linmodel operator...
            //todo: remove that fix and correct the linemodel operator to avoid this case
                if(LineModelSolutions[solutionIdx].Amplitudes[j]>0.0){
                    nSol++;
                    indexesSols.push_back(LineModelSolutions[solutionIdx].ElementId[j]);
                }
            }
            //*/
        }

    }else{
        nSol = 0;
    }

    return nSol;
}


Float64 CLineModelResult::GetExtremaMerit(Int32 extremaIdx) const
{
    Float64 outVal=-1.0;
    if(Extrema.size()>extremaIdx)
    {
        Int32 solutionIdx=0;
        for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
        {
            if(Redshifts[i2] == Extrema[extremaIdx]){
                solutionIdx = i2;
                break;
            }
        }
        outVal = ChiSquare[solutionIdx];

    }

    return outVal;
}
