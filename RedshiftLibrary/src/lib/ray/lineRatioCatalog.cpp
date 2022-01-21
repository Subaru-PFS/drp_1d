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
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/ray/lineRatioCatalog.h"

using namespace NSEpic;

CLineRatioCatalog::CLineRatioCatalog(const std::string& name, const Float64& n_sigma_support):
  m_Name(name),
  CRayCatalog(n_sigma_support)
{

}

CLineRatioCatalog::CLineRatioCatalog(const std::string& name, const CRayCatalog& lineCatalog):
  CRayCatalog(lineCatalog),
  m_Name(name)
{
  // TODO set all nominalAmplitudes to 0
  
}

/*
void CLineRatioCatalog::addLine(const std::string& name,
				const Float64& position,
				const std::string& type,
				const std::string& force,
				const std::string& profile,
				const TAsymParams& asymParams,
				const std::string& groupName,
				const Float64& nominalAmplitude,
				const std::string& velocityGroup,
				const Float64& velocityOffset,
				const bool& enableVelocityFit,
				const Int32& id)
{
  m_LineCatalog.AddRayFromParams(name,
				 position,
				 type,
				 force,
				 profile,
				 asymParams,
				 groupName,
				 nominalamplitude,
				 velocityGroup,
				 velocityOffset,
				 enableVelocityFit,
				 id);
}
*/
void CLineRatioCatalog::addVelocity(const std::string& name, const Float64& value)
{
  if(!m_Velocities.emplace(name,value).second)
    throw GlobalException(INTERNAL_ERROR,Formatter()<< "Velocity for group " << name << " already exists");
}

void CLineRatioCatalog::setPrior(const Float64& prior)
{
  m_Prior = prior;
}
void CLineRatioCatalog::setIsmIndex(const Float64& ismIndex)
{
  m_IsmIndex = ismIndex;
}
    
const Float64& CLineRatioCatalog::getVelocity(const std::string& velGroup)
{
  return m_Velocities[velGroup];
}
