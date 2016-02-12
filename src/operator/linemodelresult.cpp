#include <epic/redshift/operator/linemodelresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

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
            stream <<  Extrema[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save bic list, on 1 line
    if(bic.size()>0){
        stream <<  "#BIC for each extrema = {";
        for ( int i=0; i<bic.size(); i++)
        {
            stream <<  bic[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save posterior list, on 1 line
    if(Posterior.size()>0){
        stream <<  "#POSTERIOR for each extrema = {";
        for ( int i=0; i<Posterior.size(); i++)
        {
            stream <<  Posterior[i] << "\t";
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
}

/**
 * \brief Prints the argument store's number of redshift results in the argument stream.
 **/
Void CLineModelResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream << "LineModelResult" << "\t" << Redshifts.size() << std::endl;
}

/**
 * \brief Returns the number of strong lines in the solution that are above the given SNR thresholds.
 * Let nSol be the number of lines above a given SNR threshold. Set it as 0.
 * If the argument extremaIdx cannot be an index for Extrema, return 0.
 * Find the index for the redshift corresponding to the extremaIdx.
 * For each amplitude in the solutions:
 *   Skip if solution already checked.
 *   Skip if solution is not a strong line.
 *   If the solution error is greater than 0:
 *     Calculate the solution SNR as amplitude / noise and the fitting SNR as amplitude / fitting error.
 *     If both SNRs are above their thresholds, include the solution in the "already checked", and increase nSol by one.
 **/
Int32 CLineModelResult::GetNLinesOverCutThreshold( Int32 extremaIdx, Float64 snrThres, Float64 fitThres ) const
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
  return nSol;
}

/**
 * \brief Returns the value of the ChiSquare of the Extrema indexed by the argument extremaIdx - if it is a valid index.
 * Let result be -1.
 * If the argument extremaIdx can be an index of Extrema:
 *   Find the first Redshift corresponding to the Extrema with the index extremaIdx.
 *   Set the result to the ChiSquare with the corresponding index.
 * Return the result.
 **/
Float64 CLineModelResult::GetExtremaMerit( Int32 extremaIdx ) const
{
  Float64 outVal=-1.0;
  if( Extrema.size()>extremaIdx )
    {
      Int32 solutionIdx=-1;
      for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
        {
	  if(Redshifts[i2] == Extrema[extremaIdx]){
	    solutionIdx = i2;
	    break;
	  }
        }
      if( solutionIdx==-1 )
	{
	  return -1;
	}
      outVal = ChiSquare[solutionIdx];
    }
    return outVal;
}
