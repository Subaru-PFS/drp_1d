#include <RedshiftLibrary/operator/linemodelresult.h>

#include <RedshiftLibrary/statistics/deltaz.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <string>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/linetags.h>

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
 * @brief CLineModelResult::Init
 * Initializes the linemodel results, by:
 *  - setting the redshift grid
 *  - setting the size of chisquare vectors
 * @param redshifts
 * @param restRays
 * @param nTplshapes
 * @return
 * ERR = -2 : if the tplshapesPriors list is not the same size as the tpl-ratios
 */
Int32 CLineModelResult::Init(std::vector<Float64> redshifts,
                              CRayCatalog::TRayVector restRays,
                              Int32 nTplshapes,
                              std::vector<Float64> tplshapesPriors)
{
    Int32 err = 0;
    if(tplshapesPriors.size()!=nTplshapes)
    {
        return -2;
    }

    Int32 nResults = redshifts.size();
    Status.resize( nResults );
    ChiSquare.resize( nResults );
    ScaleMargCorrection.resize( nResults );
    Redshifts.resize( nResults );
    Redshifts = redshifts;
    restRayList = restRays;
    LineModelSolutions.resize( nResults );
    ContinuumModelSolutions.resize( nResults );

    //init the tplshape chisquare results
    for(Int32 k=0; k<nTplshapes; k++)
    {
        TFloat64List _chi2tpl(nResults, DBL_MAX);
        ChiSquareTplshapes.push_back(_chi2tpl);

        TFloat64List _corr(nResults, 0.0);
        ScaleMargCorrectionTplshapes.push_back(_corr);

        TBoolList selp(nResults, false);
        StrongELPresentTplshapes.push_back(selp);

        std::vector<Int32> nlac(nResults, false);
        NLinesAboveSNRTplshapes.push_back(nlac);

        PriorTplshapes.push_back(tplshapesPriors[k]);
        TFloat64List _logpriors(nResults, 0.0);
        PriorLinesTplshapes.push_back(_logpriors);
    }

    ChiSquareContinuum.resize( nResults );
    ScaleMargCorrectionContinuum.resize( nResults);

    return err;
}

Int32 CLineModelResult::SetChisquareTplshapeResult(Int32 index_z,
                                                    const TFloat64List & chisquareTplshape,
                                                    const TFloat64List & scaleMargCorrTplshape,
                                                    const TBoolList & strongEmissionLinePresentTplshape,
                                                    const TInt32List & nLinesAboveSNRTplshape,
                                                    const TFloat64List & priorLinesTplshape)
{
    if(index_z>=Redshifts.size())
    {
        return -1;
    }
    if(chisquareTplshape.size()!=ChiSquareTplshapes.size())
    {
        return -2;
    }
    if(chisquareTplshape.size()<1)
    {
        return -3;
    }
    if(chisquareTplshape.size()!=scaleMargCorrTplshape.size())
    {
        return -4;
    }
    if(chisquareTplshape.size()!=strongEmissionLinePresentTplshape.size())
    {
        return -4;
    }
    if(chisquareTplshape.size()!=nLinesAboveSNRTplshape.size())
    {
        return -4;
    }
    if(chisquareTplshape.size()!=priorLinesTplshape.size())
    {
        return -4;
    }

    for(Int32 k=0; k<chisquareTplshape.size(); k++)
    {
        ChiSquareTplshapes[k][index_z] = chisquareTplshape[k];
        ScaleMargCorrectionTplshapes[k][index_z] = scaleMargCorrTplshape[k];
        StrongELPresentTplshapes[k][index_z] = strongEmissionLinePresentTplshape[k];
        NLinesAboveSNRTplshapes[k][index_z] = nLinesAboveSNRTplshape[k];
        PriorLinesTplshapes[k][index_z] = priorLinesTplshape[k];
    }
    return 0;
}

TFloat64List CLineModelResult::GetChisquareTplshapeResult( Int32 index_z )
{
    TFloat64List chisquareTplshape;
    if(index_z>=Redshifts.size())
    {
        return chisquareTplshape;
    }
    if(ChiSquareTplshapes.size()<1)
    {
        return chisquareTplshape;
    }

    for(Int32 k=0; k<ChiSquareTplshapes.size(); k++)
    {
        chisquareTplshape.push_back(ChiSquareTplshapes[k][index_z]);
    }

    return chisquareTplshape;
}

TFloat64List CLineModelResult::GetScaleMargCorrTplshapeResult( Int32 index_z )
{
    TFloat64List scaleMargCorrTplshape;
    if(index_z>=Redshifts.size())
    {
        return scaleMargCorrTplshape;
    }
    if(ScaleMargCorrectionTplshapes.size()<1)
    {
        return scaleMargCorrTplshape;
    }

    for(Int32 k=0; k<ScaleMargCorrectionTplshapes.size(); k++)
    {
        scaleMargCorrTplshape.push_back(ScaleMargCorrectionTplshapes[k][index_z]);
    }

    return scaleMargCorrTplshape;
}

TBoolList CLineModelResult::GetStrongELPresentTplshapeResult( Int32 index_z )
{
    TBoolList strongELPresentTplshape;
    if(index_z>=Redshifts.size())
    {
        return strongELPresentTplshape;
    }
    if(StrongELPresentTplshapes.size()<1)
    {
        return strongELPresentTplshape;
    }

    for(Int32 k=0; k<StrongELPresentTplshapes.size(); k++)
    {
        strongELPresentTplshape.push_back(StrongELPresentTplshapes[k][index_z]);
    }

    return strongELPresentTplshape;
}


std::vector<Int32> CLineModelResult::GetNLinesAboveSNRTplshapeResult( Int32 index_z )
{
    std::vector<Int32> priorTplshape;
    if(index_z>=Redshifts.size())
    {
        return priorTplshape;
    }
    if(NLinesAboveSNRTplshapes.size()<1)
    {
        return priorTplshape;
    }

    for(Int32 k=0; k<NLinesAboveSNRTplshapes.size(); k++)
    {
        priorTplshape.push_back(NLinesAboveSNRTplshapes[k][index_z]);
    }

    return priorTplshape;
}


std::vector<Float64> CLineModelResult::GetPriorLinesTplshapeResult( Int32 index_z )
{
    std::vector<Float64> priorTplshape;
    if(index_z>=Redshifts.size())
    {
        return priorTplshape;
    }
    if(PriorLinesTplshapes.size()<1)
    {
        return priorTplshape;
    }

    for(Int32 k=0; k<PriorLinesTplshapes.size(); k++)
    {
        priorTplshape.push_back(PriorLinesTplshapes[k][index_z]);
    }

    return priorTplshape;
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
void CLineModelResult::Load( std::istream& stream )
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
            catch (boost::bad_lexical_cast&)
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
 * Print each SigmaZ as a comment.
 * Print each LogArea as a comment.
 **/
void CLineModelResult::Save( std::ostream& stream ) const
{
    stream <<  "#Redshifts\tChiSquare"<< std::endl;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        stream <<  Redshifts[i] << std::setprecision(32) << "\t" << std::scientific << ChiSquare[i] << std::fixed << std::endl;
    }

    //ExtremaResult.Save(stream); //todo: move these into their own file. aview should be adapted.

    // save dTransposeD, on 1 line
    if(true){
        stream <<  "#dTransposeD = {";
        stream <<  dTransposeD << "\t";
        stream << "}" << std::endl;
    }

}

/**
 * \brief Prints the argument store's number of redshift results in the argument stream.
 **/
void CLineModelResult::SaveLine(  std::ostream& stream ) const
{
    stream << "LineModelResult" << "\t" << Redshifts.size() << std::endl;
}

Int32 CLineModelResult::GetNLinesOverCutThreshold(Int32 solutionIdx, Float64 snrThres, Float64 fitThres) const
{
    Int32 nSol=0;

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
 * @brief CLineModelResult::GetStrongLinesPresence
 * @param filterType: 1: emission only, 2 abs only, else: no filter
 * @return: a list of boolean values indicating if a strong is present (not outsidelambdarange for that z) for each redshift
 */
TBoolList CLineModelResult::GetStrongLinesPresence( UInt32 filterType, std::vector<CLineModelSolution> linemodelsols ) const
{
    TBoolList strongIsPresent(linemodelsols.size(), false);
    for ( UInt32 solutionIdx=0; solutionIdx<linemodelsols.size(); solutionIdx++)
    {
        strongIsPresent[solutionIdx] = false;

        for ( UInt32 j=0; j<linemodelsols[solutionIdx].Amplitudes.size(); j++)
        {
            if( !linemodelsols[solutionIdx].Rays[j].GetIsStrong() )
            {
                continue;
            }

            if(filterType==1)
            {
                if( !linemodelsols[solutionIdx].Rays[j].GetIsEmission() )
                {
                    continue;
                }
            }else if(filterType==2)
            {
                if( linemodelsols[solutionIdx].Rays[j].GetIsEmission() )
                {
                    continue;
                }
            }

            if( linemodelsols[solutionIdx].OutsideLambdaRange[j] )
            {
                continue;
            }

            if(linemodelsols[solutionIdx].Amplitudes[j]>0.0)
            {
                strongIsPresent[solutionIdx] = true;
                break;
            }
        }
    }

    return strongIsPresent;
}


std::vector<Int32> CLineModelResult::GetNLinesAboveSnrcut( std::vector<CLineModelSolution> linemodelsols ) const
{
    std::vector<Int32> nlinesabove(linemodelsols.size(), 0);
    for ( UInt32 solutionIdx=0; solutionIdx<linemodelsols.size(); solutionIdx++)
    {
        nlinesabove[solutionIdx] = linemodelsols[solutionIdx].NLinesAboveSnrCut;
    }


    return nlinesabove;
}


/**
 * WARNING: this function has not been tested at all !!! please check/debug
 * @brief CLineModelResult::GetStrongestLineIsHa
 * @return: a list of boolean values indicating if the strongest line is Ha (Highest amp and not outsidelambdarange for that z) for each redshift
 */
TBoolList CLineModelResult::GetStrongestLineIsHa( std::vector<CLineModelSolution> linemodelsols ) const
{
    Bool verbose = true;
    UInt32 filterType=1;
    linetags ltags;
    TBoolList strongestIsHa(linemodelsols.size(), false);
    std::string ampMaxLineTag = "";
    for ( UInt32 solutionIdx=0; solutionIdx<linemodelsols.size(); solutionIdx++)
    {
        strongestIsHa[solutionIdx] = false;
        Float64 ampMax = -1;
        Float64 ampHa = -1;
        for ( UInt32 j=0; j<linemodelsols[solutionIdx].Amplitudes.size(); j++)
        {

            if(filterType==1)
            {
                if( !linemodelsols[solutionIdx].Rays[j].GetIsEmission() )
                {
                    continue;
                }
            }else if(filterType==2)
            {
                if( linemodelsols[solutionIdx].Rays[j].GetIsEmission() )
                {
                    continue;
                }
            }

            if( linemodelsols[solutionIdx].OutsideLambdaRange[j] )
            {
                continue;
            }

            Log.LogDebug("    linemodelresult: using ray for max amp search=%s", linemodelsols[solutionIdx].Rays[j].GetName().c_str());
            if(linemodelsols[solutionIdx].Amplitudes[j]>ampMax)
            {
                ampMax = linemodelsols[solutionIdx].Amplitudes[j];
                ampMaxLineTag = linemodelsols[solutionIdx].Rays[j].GetName().c_str();
            }
            if(linemodelsols[solutionIdx].Rays[j].GetName()==ltags.halpha_em)
            {
                ampHa = linemodelsols[solutionIdx].Amplitudes[j];
            }
        }

        if(ampHa>0 && ampMax==ampHa)
        {
            strongestIsHa[solutionIdx] = true;
        }
        if(verbose)
        {
            Log.LogDetail("    linemodelresult: z=%f, ampHa=%e, ampMax=%e, ampMaxLineTag=%s", linemodelsols[solutionIdx].Redshift, ampHa, ampMax, ampMaxLineTag.c_str());
        }
    }

    return strongestIsHa;
}

Float64 CLineModelResult::GetMinChiSquare() const
{
    Float64 min=DBL_MAX;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        if(min>ChiSquare[i])
        {
            min= ChiSquare[i];
        }
    }
    return min;
}

Float64 CLineModelResult::GetMaxChiSquare() const
{
    Float64 max=-DBL_MAX;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        if(max<ChiSquare[i])
        {
            max= ChiSquare[i];
        }
    }
    return max;
}

Int32 CLineModelResult::getRedshiftIndex(Float64 z)
{
  std::vector<Float64>::iterator itr = std::lower_bound(Redshifts.begin(),
                                                        Redshifts.end(),
                                                        z);

  if (itr == Redshifts.end() || *itr != z)
    {
      Log.LogError("CLineModelResult::getRedshiftIndex: Could not find redshift index for %f", z);
      throw runtime_error("CLineModelResult::getRedshiftIndex: Could not find redshift index");
    }
  return (itr - Redshifts.begin()); 
}
