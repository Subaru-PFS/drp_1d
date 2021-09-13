// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/linemodel/element.h"

#include "RedshiftLibrary/debug/assert.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/exception.h"


#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <vector>

namespace bfs = boost::filesystem;


using namespace NSEpic;

CLineModelElement::CLineModelElement(const std::string& widthType, const Float64 velocityEmission, const Float64 velocityAbsorption):
    m_VelocityEmission(velocityEmission),
    m_VelocityAbsorption(velocityAbsorption),
    m_dataExtinctionFlux(NULL),
    m_OutsideLambdaRange(true),
    m_fittingGroupInfo("-1"),
    m_OutsideLambdaRangeOverlapThreshold(0.33)//33% overlap minimum in order to keep the line
    //example: 0.33 means 66% of the line is allowed to be outside the spectrum with the line still considered inside the lambda range

{
    if( widthType == "instrumentdriven"){
        m_LineWidthType = INSTRUMENTDRIVEN;
    } else if( widthType == "combined"){
        m_LineWidthType = COMBINED;
    } else if( widthType == "velocitydriven"){
        m_LineWidthType = VELOCITYDRIVEN;
    }else {
        Log.LogError("Unknown LineWidthType %s", widthType.c_str());
        throw std::runtime_error("Unknown LineWidthType");
    }

    //LoadDataExtinction(); //uncomment if this line profile is used
}

CLineModelElement::~CLineModelElement()
{
}

std::string CLineModelElement::GetElementTypeTag()
{
    return m_ElementType;
}

Int32 CLineModelElement::FindElementIndex(Int32 LineCatalogIndex)
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_LineCatalogIndexes.size(); iElts++ )
    {
        if(m_LineCatalogIndexes[iElts] == LineCatalogIndex){
            idx = iElts;
            break;
        }
    }

    return idx;
}

Int32 CLineModelElement::GetSize()
{
    return (Int32)m_LineCatalogIndexes.size();
}

bool CLineModelElement::IsOutsideLambdaRange()
{
    return m_OutsideLambdaRange;
}
//redirecting to the lsf method for computing instrument responce
/**
 * Get instrumental response (including source response) from LSF
 * combine quadratically with the instrinsic width of the Line itself. The line width in this case represents to the velocity
 * */
Float64 CLineModelElement::GetLineWidth(Float64 redshiftedlambda, Float64 z, Bool isEmission)
{
    const Float64 c = m_speedOfLightInVacuum;
    Float64 v = isEmission ? m_VelocityEmission : m_VelocityAbsorption;
    const Float64 pfsSimuCompensationFactor = 1.0;

    if(!m_LSF)
        throw std::runtime_error("LSF object is not initailized.");
    Float64 instrumentSigma = m_LSF->GetWidth(redshiftedlambda);
   
    Float64 velocitySigma = pfsSimuCompensationFactor*v/c*redshiftedlambda;//, useless /(1+z)*(1+z);
    switch(m_LineWidthType)
    {
        case INSTRUMENTDRIVEN://only instrumental sigma
            velocitySigma=0.;
            break;
        case VELOCITYDRIVEN://only velocity sigma
            instrumentSigma = 0.;
            break;
        case COMBINED://combination of the two
            break;
        default:
            Log.LogError("Invalid LSFType %d", m_LineWidthType);
            throw std::runtime_error("Unknown mode");
    }

    Float64 sigma = sqrt(instrumentSigma*instrumentSigma + velocitySigma*velocitySigma);
    return sigma;
}
Float64 CLineModelElement::GetLineProfileDerivVel(std::shared_ptr<CLineProfile>& profile, Float64 x, Float64 x0, Float64 sigma, Bool isEmission)
{
    const Float64 c = m_speedOfLightInVacuum;
    const Float64 pfsSimuCompensationFactor = 1.0;
    Float64 v = isEmission ? m_VelocityEmission : m_VelocityAbsorption,
            v_to_sigma = pfsSimuCompensationFactor/c*x0;
    
    //sincs lsf is an instrumental response, then derivative of this latter with respect to velocity is null   
    switch (m_LineWidthType) {
        case INSTRUMENTDRIVEN:
        //case FIXED:
            return 0.0;
        case COMBINED:
    //    case NISPSIM2016:
    //    case NISPVSSPSF201707: //not supported as of 2017-07
            return v_to_sigma * v_to_sigma * v /sigma * profile->GetLineProfileDerivSigma( x, x0, sigma);
        case VELOCITYDRIVEN:
            return v_to_sigma * profile->GetLineProfileDerivSigma( x, x0, sigma);
        default:
            Log.LogError("Invalid LineWidthType : %d", m_LineWidthType);
            throw std::runtime_error("Unknown LineWidthType");
    }
    return 0.0;
}
//TODO: check if below can be removed
//this is called from CLineModelElementList: thus the call should be redirected to the new lsf class
void CLineModelElement::SetSourcesizeDispersion(Float64 sigma)
{
    //m_SourceSizeDispersion = sigma;
    m_LSF->SetSourcesizeDispersion(sigma);
}

void CLineModelElement::SetLSF(const std::shared_ptr<const CLSF> & lsf)
{
    m_LSF = lsf;
}

void CLineModelElement::SetVelocityEmission(Float64 vel)
{
    m_VelocityEmission = vel;
}

void CLineModelElement::SetVelocityAbsorption(Float64 vel)
{
    m_VelocityAbsorption = vel;
}

Float64 CLineModelElement::GetVelocityEmission()
{
    return m_VelocityEmission;
}

Float64 CLineModelElement::GetVelocityAbsorption()
{
    return m_VelocityAbsorption;
}

Float64 CLineModelElement::GetVelocity()
{
    Float64 vel=-1;
    if(m_Rays.size()>0)
    {
        if(m_Rays[0].GetIsEmission())
        {
            vel = m_VelocityEmission;
        }else{
            vel = m_VelocityAbsorption;
        }
    }
    return vel;
}

void CLineModelElement::setVelocity(Float64 vel)
{
    if(m_Rays.size()>0)
    {
        if(m_Rays[0].GetIsEmission())
        {
            m_VelocityEmission = vel;
        }else{
            m_VelocityAbsorption = vel;
        }
    }
    else throw GlobalException(INTERNAL_ERROR,"Empty line model element, could not set velocity");

}

//wrapper function 
void CLineModelElement::SetAsymfitParams(TAsymParams params, Int32 idx)
{
    if(!m_asymLineIndices.size()) return;
    if(idx>=0){
        m_Rays[idx].SetAsymParams(params);
    }
    else{
        for(auto i : m_asymLineIndices)
            m_Rays[i].SetAsymParams(params);
    }
    return;
}
//wrapper function 
void CLineModelElement::resetAsymfitParams()
{
    for(auto i : m_asymLineIndices)
        m_Rays[i].resetAsymFitParams();
}

const TAsymParams CLineModelElement::GetAsymfitParams(UInt32 idx)
{
    if(!m_asymLineIndices.size())
        return {NAN, NAN, NAN};//case where no asymprofile in linecatalog
    return m_Rays[m_asymLineIndices[idx]].GetAsymParams();
}

Float64 CLineModelElement::GetSumCross()
{
    return m_sumCross;
}

void CLineModelElement::SetSumCross(Float64 val)
{
    m_sumCross=val;
}

Float64 CLineModelElement::GetDtmFree()
{
    return m_dtmFree;
}

void CLineModelElement::SetDtmFree(Float64 val)
{
    m_dtmFree=val;
}

Float64 CLineModelElement::GetSumGauss()
{
    return m_sumGauss;
}

void CLineModelElement::SetSumGauss(Float64 val)
{
    m_sumGauss=val;
}

Float64 CLineModelElement::GetFitAmplitude()
{
    return m_fitAmplitude;
}


/**
 * Load the extinction residue data.
 */
Bool CLineModelElement::LoadDataExtinction()
{
    std::string filePathStr = "/home/aschmitt/Documents/amazed/methods/linemodel/extinction_element/extinction_data/element_meiksin_dl0.1cropped930-1230.fits.txt";
    Log.LogDebug ( "Parsing ASCII file %s.", filePathStr.c_str() );
    if( !bfs::exists( filePathStr.c_str() ) )
    {
        Log.LogError( "Read: Path for element_extinction file does not exist." );
        return false;
    }

    bfs::ifstream file;
    file.open( filePathStr.c_str() );

    m_dataExtinctionFlux = (Float64 *) malloc(m_dataN* sizeof(Float64));

    Int32 i = 0;
    file.clear();
    file.seekg( 0 );

    for( std::string line; std::getline( file, line ); )
    {
        if( !boost::starts_with( line, "#" ) )
        {
            std::istringstream iss( line );
            Float64 x, y;
            iss >> x >> y;
            //spcSpectralAxis[i] = x; //not reading x, hardcoded from 930 to 1230, with dlambda = 0.1A
            m_dataExtinctionFlux[i] = y;
            i++;
            if(i>=m_dataN)
            {
                break;
            }
        }
    }
    file.close();
    if(i<m_dataN)
    {
        Log.LogError( "Read: linemodel- extinction residue: data not read successfully" );
    }
    Log.LogDetail( "Read exctinction successfully" );
    return true;
}
