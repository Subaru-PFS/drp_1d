#include "RedshiftLibrary/linemodel/modelfittingresult.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"

using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CModelFittingResult::CModelFittingResult()
{
  this->m_type = "CModelFittingResult";
}


/**
 * \brief Attributes values to member variables according to arguments.
 **/
CModelFittingResult::CModelFittingResult( CLineModelSolution _lineModelSolution,
                                          Float64 _redshift,
                                          Float64 _merit,
                                          CRayCatalog::TRayVector _restRayList,
                                          Float64 _velEmission,
                                          Float64 _velAbsorption)
{
    this->m_type = "CModelFittingResult";
    LineModelSolution   = _lineModelSolution;
    Redshift            = _redshift;
    Merit               = _merit;
    restRayList         = _restRayList;
    VelocityEmission               = _velEmission;
    VelocityAbsorption               = _velAbsorption;

    // This is only to make linemodelsolution members available for output
    for (UInt32 j=0; j<LineModelSolution.Amplitudes.size(); j++) rayId.emplace_back(restRayList[j].GetID());
    for (UInt32 j=0; j<LineModelSolution.Amplitudes.size(); j++) rayLambdaRest.emplace_back(restRayList[j].GetPosition());
    FittedRaysFlux = LineModelSolution.Fluxs;
    FittedRaysLambda = LineModelSolution.LambdaObs;    
}

/**
 * \brief Empty destructor.
 **/
CModelFittingResult::~CModelFittingResult()
{
}


//Load the linemodel fit results from a csv file;
//WARNING: read only the amplitudes fitted so far, as of 2016-03-10
void CModelFittingResult::Load( const char* filePath )
{
    std::ifstream file;
    file.open( filePath, std::ifstream::in );
    if( file.rdstate() & std::ios_base::failbit )
    {
        return;
    }

    std::string line;

    // Read file line by line
    while( getline( file, line ) )
    {
        if( boost::starts_with( line, "#" ) )
        {
            continue;
        }
        boost::char_separator<char> sep("\t");
        // Tokenize each line
        typedef boost::tokenizer< boost::char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );

        // Check if it's not a comment
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() )
        {
            // skip type
            ++it;
            // skip force
            ++it;

            //parse name
            //std::string name;
            if( it != tok.end() )
            {
              // name = *it;
            }
            else
            {
	      file.close();
	      return;
            }

            // skip
            ++it;
            // skip
            ++it;

            ++it;
            // Parse amplitude fitted
            Float64 amp = 0.0;
            try
            {
                amp = boost::lexical_cast<double>(*it);
            }
            catch (boost::bad_lexical_cast&)
            {
                amp = 0.0;
		file.close();
		return;
            }
            LineModelSolution.Amplitudes.push_back(amp);
	    

        }
    }
    file.close();
}
