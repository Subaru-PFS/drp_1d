#include <RedshiftLibrary/linemodel/modelcontinuumfittingresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <RedshiftLibrary/spectrum/spectrum.h>

using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CModelContinuumFittingResult::CModelContinuumFittingResult()
{
}


/**
 * \brief Attributes values to member variables according to arguments.
 **/
CModelContinuumFittingResult::CModelContinuumFittingResult( Float64 _redshift,
                                                            std::string _name,
                                                            Float64 _merit,
                                                            Float64 _amp,
                                                            Float64 _ismCoeff,
                                                            Int32 _igmIndex)
{
    Redshift            = _redshift;
    Merit               = _merit;
    Amp               = _amp;
    Name               = _name;
    IsmCoeff               = _ismCoeff;
    IgmIndex               = _igmIndex;

}

/**
 * \brief Empty destructor.
 **/
CModelContinuumFittingResult::~CModelContinuumFittingResult()
{
}

/**
 * \brief Prints the results of the Linemodel in the argument store, using the argument stream as output.
 **/
void CModelContinuumFittingResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    // save model continuum solution

    stream <<  "#z\ttpl_name\tmerit\tamp\tism_coeff\tigm_index" <<  std::endl;
    stream << std::fixed <<  Redshift << "\t";
    stream << Name.c_str() << "\t";
    stream << std::scientific << std::setprecision(5) << Merit << "\t";
    stream << std::scientific << std::setprecision(5) << Amp << "\t";
    stream << std::fixed << std::setprecision(3) << IsmCoeff << "\t";
    stream << std::fixed << std::setprecision(1) << IgmIndex <<  std::endl;

}

/**
 * \brief Empty method.
 **/
void CModelContinuumFittingResult::SaveLine(const CDataStore &store, std::ostream& stream ) const
{

}
