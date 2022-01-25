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
#ifndef _LINE_RATIO_CATALOG_H_
#define _LINE_RATIO_CATALOG_H_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/ray.h"
#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

#include <boost/format.hpp>

#include <vector>
#include <string>

namespace NSEpic
{

  class CLineRatioCatalog : public CRayCatalog
  {
  public:
    CLineRatioCatalog(const std::string& name, const Float64& n_sigma_support);
    /*CLineRatioCatalog(const std::string& name,
		      const Float64& n_sigma_support,
		      const Float64& prior
		      const Int32& ismIndex);*/
    CLineRatioCatalog(const std::string& name, const CRayCatalog& lineCatalog);
    ~CLineRatioCatalog() = default;
    
    /*    void addLine(const std::string& name,
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
		 const Int32& id);
    */

    void addVelocity(const std::string& name, const Float64& value);
    void setPrior(const Float64& prior);
    void setIsmIndex(const Float64& ismIndex);

    const Float64& getPrior(){return m_Prior;}
    const std::string& getName(){return m_Name;}
    const Int32& getIsmIndex(){ return m_IsmIndex;}
    const Float64& getVelocity(const std::string&group);
  private:
    //    CRayCatalog m_LineCatalog;
    std::string m_Name;
    std::map<std::string, Float64> m_Velocities;
    Float64 m_Prior;
    Int32 m_IsmIndex = 0;
    
    
  };
}

#endif
