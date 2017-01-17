#include <epic/redshift/linemodel/modelfittingresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <epic/redshift/spectrum/spectrum.h>

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
CModelFittingResult::CModelFittingResult( CLineModelResult::SLineModelSolution _lineModelSolution,
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
Void CModelFittingResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    // save linemodel solution
    if(true){
        for ( UInt32 i=0; i<1; i++)
        {
            stream <<  "#linemodel solution " << i << " for z = " <<  std::fixed <<  Redshift;
            stream << ", velocityEmission = " <<  VelocityEmission;
            stream << ", velocityAbsorption = " <<  VelocityAbsorption;
            stream << ", merit = " <<  Merit << "{" <<  std::endl;
            stream << "#type\t#force\t#Name_____________\t#elt_ID\t#lambda_rest\t#amp_____\t#err_____\t#err_fit_____\n";
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
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.Amplitudes[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.Errors[j] << "\t";
                stream << std::scientific << std::setprecision(5) <<  LineModelSolution.FittingError[j] << std::endl;
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
Void CModelFittingResult::Load( const char* filePath )
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
            std::string name;
            if( it != tok.end() )
            {
                name = *it;
            }
            else
            {
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
            catch (boost::bad_lexical_cast)
            {
                amp = 0.0;
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
Void CModelFittingResult::SaveLine(const CDataStore &store, std::ostream& stream ) const
{

}
