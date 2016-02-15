#include <epic/redshift/linemodel/modelfittingresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
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
CModelFittingResult::CModelFittingResult(CLineModelResult::SLineModelSolution _lineModelSolution, Float64 _redshift, Float64 _merit, CRayCatalog::TRayVector _restRayList)
{
    LineModelSolution   = _lineModelSolution;
    Redshift            = _redshift;
    Merit               = _merit;
    restRayList         = _restRayList;
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
        }
    }

}

/**
 * \brief Empty method.
 **/
Void CModelFittingResult::SaveLine(const CDataStore &store, std::ostream& stream ) const
{
}
