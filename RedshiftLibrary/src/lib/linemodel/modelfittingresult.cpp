#include <RedshiftLibrary/linemodel/modelfittingresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/log/log.h>

using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CModelFittingResult::CModelFittingResult()
{
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
    LineModelSolution   = _lineModelSolution;
    Redshift            = _redshift;
    Merit               = _merit;
    restRayList         = _restRayList;
    VelocityEmission               = _velEmission;
    VelocityAbsorption               = _velAbsorption;

}

/**
 * \brief Empty destructor.
 **/
CModelFittingResult::~CModelFittingResult()
{
}

/**
 * \brief Prints the results of the Linemodel in the argument store, using the argument stream as output.
 **/
void CModelFittingResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    // save linemodel solution
    if(true){
        for ( UInt32 i=0; i<1; i++)
        {
            stream <<  "#linemodel solution " << i << " for z = " <<  std::fixed <<  Redshift;
            stream << ", velocityEmission = " <<  VelocityEmission;
            stream << ", velocityAbsorption = " <<  VelocityAbsorption;
            stream << ", merit = " <<  Merit << "{" <<  std::endl;
            stream << "#type\t#force\t#Name_____________\t#elt_ID\t#lambda_rest_beforeOffset\t#lambda_obs\t#amp_____\t#err_____\t#err_fit_____\t#fit_group_____\t#velocity_____\t#offset_____\t#sigma_____\t#flux_____\t#flux_err_____\t#flux_di_____\t#center_cont_flux_____\t#cont_err_____\n";
            for ( UInt32 j=0; j<LineModelSolution.Amplitudes.size(); j++)
            {
                std::string typeStr="";
                if(restRayList[j].GetType() == CRay::nType_Absorption){
                    typeStr = "A";
                }else{
                    typeStr = "E";
                }
                stream <<  typeStr << "\t";
                std::string forceStr="";
                if(restRayList[j].GetForce() == CRay::nForce_Strong){
                    forceStr = "S";
                }else{
                    forceStr = "W";
                }
                stream <<  forceStr << "\t";
                std::string name = restRayList[j].GetName();
                Int32 nstr = name.size();
                for(int jstr=0; jstr<18-nstr; jstr++){
                    name = name.append(" ");
                }
                stream <<  std::fixed << name << "\t";
                stream <<  std::fixed << std::setprecision(0) << LineModelSolution.ElementId[j] << "\t";
                stream <<  std::fixed << std::setprecision(3) << restRayList[j].GetPosition() << "\t";
                stream <<  std::fixed << std::setprecision(3) << LineModelSolution.LambdaObs[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.Amplitudes[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.Errors[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.FittingError[j] << "\t";
                stream << std::fixed << std::setprecision(1) <<  LineModelSolution.fittingGroupInfo[j] << "\t";
                stream << std::fixed << std::setprecision(1) <<  LineModelSolution.Velocity[j] << "\t";
                stream << std::fixed << std::setprecision(1) <<  LineModelSolution.Offset[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.Sigmas[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.Fluxs[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.FluxErrors[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.FluxDirectIntegration[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.CenterContinuumFlux[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.ContinuumError[j] << std::endl;
            }
            stream << "#}" << std::endl;

            stream <<  "#lya asymfit params = { ";
            stream << "widthCoeff= " <<  LineModelSolution.LyaWidthCoeff << ", ";
            stream << "alpha= " <<  LineModelSolution.LyaAlpha << ", ";
            stream << "delta= " <<  LineModelSolution.LyaDelta << " ";
            stream << "}" << std::endl;
        }
    }

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


/**
 * \brief Empty method.
 **/
void CModelFittingResult::SaveLine(const CDataStore &store, std::ostream& stream ) const
{

}


void CModelFittingResult::getData(const std::string& name, int **data, int *size) const
{
  *size = LineModelSolution.LambdaObs.size();
  //  std::vector<int> rayId = std::vector<int>(LineModelSolution.LambdaObs.size()); // beware of memory leaks here, who destroys this vector ?
 
  if (name.compare("FittedRaysID") == 0)
    {

            if (rayId.empty())
      	{
          //      	  rayId = std::vector<int>(*size);
	  for (UInt32 j=0; j<LineModelSolution.Amplitudes.size(); j++) rayId.emplace_back(restRayList[j].GetID());
        }
	  //	}
	  *data = const_cast<int *>(rayId.data());
    }
    else Log.LogError("unkwown data %s",name.c_str());

}


void CModelFittingResult::getData(const std::string& name, double **data, int *size) const
{
  *size = LineModelSolution.LambdaObs.size();
  if (name.compare("FittedRaysLambda") == 0)
    {
      *data = const_cast<double *>(LineModelSolution.LambdaObs.data());
    }
  /*  else if (name.compare("FittedRaysLambdaRest") == 0)
    {
      if (lambdaRest.empty())
	{
	  lambdaRest = new std::vector<double>(*size);
	  for (UInt32 j=0; j<LineModelSolution.Amplitudes.size(); j++) lambdaRest.emplace_back(restRayList[j].GetPosition());
	}
      *data = const_cast<double *>(lambdaRest.data());
      }*/
  else if (name.compare("FittedRaysFlux") == 0)
    {
      *data = const_cast<double *>(LineModelSolution.Fluxs.data());
    }
  else Log.LogError("unkwown data %s",name.c_str());

}

void CModelFittingResult::getData(const std::string& name, std::string *data, int *size) const
{
  *size  = LineModelSolution.LambdaObs.size();
  /*
 if (name.compare("FittedRaysType") == 0)
    {
      if (rayType.empty())
	{
	  rayType = new std::vector<std::string>(*size);
	  for (UInt32 j=0; j<LineModelSolution.Amplitudes.size(); j++)
	    {
	      if(restRayList[j].GetType() == CRay::nType_Absorption) rayType.emplace_back("A");
	      else rayType.emplace_back("E");
	    }
	}
      *data = const_cast<double *>(rayType.data());
    }
  else if (name.compare("FittedRaysForce") == 0)
    {
      if (rayForce.empty())
	{
	  rayForce = new std::vector<std::string>(*size);
	  for (UInt32 j=0; j<LineModelSolution.Amplitudes.size(); j++)
	    {
	      if(restRayList[j].GetForce() == CRay::nForce_Strong) rayType.emplace_back("S");
	      else rayType.emplace_back("W");
	    }
	}
      *data = const_cast<double *>(rayType.data());
    }
  else if (name.compare("FittedRaysName") == 0)
    {
      if (rayForce.empty())
	{
	  rayForce = new std::vector<std::string>(*size);
	  for (UInt32 j=0; j<LineModelSolution.Amplitudes.size(); j++)
	    {
	      if(restRayList[j].GetForce() == CRay::nForce_Strong) rayType.emplace_back("S");
	      else rayType.emplace_back("W");
	    }
	}
      *data = const_cast<double *>(rayType.data());
    }
  */ 
  

}
