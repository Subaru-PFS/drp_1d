#include <epic/redshift/operator/linemodelresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <string>

using namespace NSEpic;


/**
 * \brief Empty constructor.
 **/
CLineModelResult::CLineModelResult()
{

}

/**
 * \brief Empty destructor.
 **/
CLineModelResult::~CLineModelResult()
{

}

Void CLineModelResult::ResizeExtremaResults(Int32 size)
{
    Extrema.resize(size);
    ExtremaMerit.resize(size);
    Posterior.resize(size);
    StrongELSNR.resize(size);
    LogArea.resize(size);
    LogAreaCorrectedExtrema.resize(size);
    SigmaZ.resize(size);
    bic.resize(size);
    ContinuumIndexes.resize(size);
    OutsideLinesMask.resize(size);
    FittedTplName.resize(size);
    FittedTplAmplitude.resize(size);
    FittedTplcorrTplName.resize(size);
}


/**
 * \brief Attempts to read the redshift and chisquare values stored in the argument stream.
 * Clear the current lines lists.
 * For each line in input stream:
 *   Tokenize each line.
 *   If the first character is not "#":
 *     Try to read the redshift. If this fails, return.
 *     Ignore the next token - it should be a name (string).
 *     Try to read the chisquare. If this fails, return.
 **/
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

	    // Parse c
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

/**
 * \brief Prints the results currently in the argument store, in the argument stream.
 * Using the argument stream to print values from the argument store:
 * Print a header as a comment.
 * Print each redshift and chisquare values.
 * Print each Extrema as a comment.
 * Print each BIC as a comment.
 * Print each POSTERIOR as a comment.
 * Print each SigmaZ as a comment.
 * Print each LogArea as a comment.
 **/
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

    // save ContinuumIndexes list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#ContinuumIndexes Color for each extrema = {";
        for ( int i=0; i<ContinuumIndexes.size(); i++)
        {
            stream << "<";
            for(Int32 kci=0; kci<ContinuumIndexes[i].size(); kci++)
            {
                stream <<  ContinuumIndexes[i][kci].Color << "\t";
            }
            stream << ">";
        }
        stream << "}" << std::endl;
        stream <<  "#ContinuumIndexes Break for each extrema = {";
        for ( int i=0; i<ContinuumIndexes.size(); i++)
        {
            stream << "<";
            for(Int32 kci=0; kci<ContinuumIndexes[i].size(); kci++)
            {
                stream <<  ContinuumIndexes[i][kci].Break << "\t";
            }
            stream << ">";
        }
        stream << "}" << std::endl;
    }


    // save StrongELSNR list, on 1 line
    if(StrongELSNR.size()>0){
        stream <<  "#StrongELSNR for each extrema = {";
        for ( int i=0; i<StrongELSNR.size(); i++)
        {
            stream <<  StrongELSNR[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save dTransposeDNocontinuum, on 1 line
    if(StrongELSNR.size()>0){
        stream <<  "#dTransposeDNocontinuum for each extrema = {";
        for ( int i=0; i<StrongELSNR.size(); i++)
        {
            stream <<  dTransposeDNocontinuum << "\t";
        }
        stream << "}" << std::endl;
    }

    // save dTransposeD, on 1 line
    if(StrongELSNR.size()>0){
        stream <<  "#dTransposeD for each extrema = {";
        for ( int i=0; i<StrongELSNR.size(); i++)
        {
            stream <<  dTransposeD << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplName, on 1 line
    if(FittedTplName.size()>0){
        stream <<  "#FittedTplName for each extrema = {";
        for ( int i=0; i<FittedTplName.size(); i++)
        {
            stream <<  FittedTplName[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplAmplitude, on 1 line
    if(FittedTplAmplitude.size()>0){
        stream <<  "#FittedTplAmplitude for each extrema = {";
        for ( int i=0; i<FittedTplAmplitude.size(); i++)
        {
            stream <<  FittedTplAmplitude[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplcorrTplName, on 1 line
    if(FittedTplcorrTplName.size()>0){
        stream <<  "#FittedTplcorrTplName for each extrema = {";
        for ( int i=0; i<FittedTplcorrTplName.size(); i++)
        {
            stream <<  FittedTplcorrTplName[i] << "\t";
        }
        stream << "}" << std::endl;
    }

}

/**
 * \brief Prints the argument store's number of redshift results in the argument stream.
 **/
Void CLineModelResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream << "LineModelResult" << "\t" << Redshifts.size() << std::endl;
}

Int32 CLineModelResult::GetNLinesOverCutThreshold(Int32 extremaIdx, Float64 snrThres, Float64 fitThres) const
{
    if( Extrema.size()<=extremaIdx )
    {
        return 0;
    }
    Int32 nSol=0;
    Int32 solutionIdx=0;
    for ( UInt32 i2=0; i2<LineModelSolutions.size(); i2++)
    {
        if( Redshifts[i2]==Extrema[extremaIdx] )
        {
            solutionIdx = i2;
            break;
        }
    }
    std::vector<Int32> indexesSols;
    for ( UInt32 j=0; j<LineModelSolutions[solutionIdx].Amplitudes.size(); j++)
    {
        //skip if already sol
        bool alreadysol = false;
        for( Int32 i=0; i<indexesSols.size(); i++ )
        {
            if( LineModelSolutions[solutionIdx].ElementId[j]==indexesSols[i] )
            {
                alreadysol=true;
                break;
            }
        }
        if( alreadysol )
        {
            continue;
        }
        if( !LineModelSolutions[solutionIdx].Rays[j].GetIsStrong() )
        {
            continue;
        }
        if( !LineModelSolutions[solutionIdx].Rays[j].GetIsEmission() )
        {
            continue;
        }

        Float64 noise = LineModelSolutions[solutionIdx].Errors[j];
        if( noise>0 )
        {
            Float64 snr = LineModelSolutions[solutionIdx].Amplitudes[j]/noise;
            Float64 Fittingsnr = LineModelSolutions[solutionIdx].Amplitudes[j]/LineModelSolutions[solutionIdx].FittingError[j];
            if( snr>=snrThres && Fittingsnr>=fitThres )
            {
                nSol++;
                indexesSols.push_back(LineModelSolutions[solutionIdx].ElementId[j]);
            }
        }

    }
    return nSol;
}

/**
 * \brief Returns the value of the ChiSquare of the Extrema indexed by the argument extremaIdx - if it is a valid index.
 * Let result be -1.
 * Return the result.
 **/
Float64 CLineModelResult::GetExtremaMerit( Int32 extremaIdx ) const
{
    Float64 outVal=-1.0;
    if( Extrema.size()>extremaIdx && ExtremaMerit.size()>extremaIdx )
    {
        outVal = ExtremaMerit[extremaIdx];
    }
    return outVal;
}
