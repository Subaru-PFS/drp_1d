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
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/ray/catalog.h"

#include <algorithm>    // std::sort
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <boost/filesystem.hpp>

#include "RedshiftLibrary/ray/lineprofile.h"

using namespace NSEpic;
using namespace std;
using namespace boost;
using namespace boost::filesystem;


CRayCatalog::CRayCatalog():
  m_nSigmaSupport(N_SIGMA_SUPPORT)
{

}

CRayCatalog::CRayCatalog(Float64 sigmaSupport):
  m_nSigmaSupport(sigmaSupport)
{

}




CRayCatalog::~CRayCatalog()
{

}

const CRayCatalog::TRayVector& CRayCatalog::GetList() const
{
    return m_List;
}

const CRayCatalog::TRayVector CRayCatalog::GetFilteredList(Int32 typeFilter, Int32 forceFilter) const
{
    {
        TRayVector filteredList;
        for( int i = 0; i< m_List.size(); i++ )
        {
            if( typeFilter == -1 || typeFilter == m_List[i].GetType()){
                if( forceFilter == -1 || forceFilter == m_List[i].GetForce()){
                    filteredList.push_back(m_List[i]);
                }
            }
        }
        return filteredList;
    }
}

const std::vector<CRayCatalog::TRayVector> CRayCatalog::ConvertToGroupList( TRayVector filteredList ) const
{

    std::vector<std::string> tags;
    for( int i = 0; i< filteredList.size(); i++ )
    {
        if(filteredList[i].GetGroupName() != "-1"){
            tags.push_back(filteredList[i].GetGroupName());
        }
    }

    // create the group tag set by removing duplicates
    std::sort( tags.begin(), tags.end() );
    tags.erase( std::unique( tags.begin(), tags.end() ), tags.end() );

    //get all group tags
    std::vector<TRayVector> fullList;
    for( Int32 itag = 0; itag<tags.size(); itag++){
        TRayVector taggedGroupList;
        for( int i = 0; i< filteredList.size(); i++ )
        {
            std::string group = filteredList[i].GetGroupName();
            if(group==tags[itag]){
                taggedGroupList.push_back(filteredList[i]);
            }
        }
        fullList.push_back(taggedGroupList);
    }
    //add the non grouped lines
    for( int i = 0; i< filteredList.size(); i++ )
    {
        std::string group = filteredList[i].GetGroupName();
        if(group=="-1"){
            TRayVector taggedGroupList;
            taggedGroupList.push_back(filteredList[i]);
            fullList.push_back(taggedGroupList);
        }
    }

    return fullList;
}


void CRayCatalog::Add( const CRay& r )
{
    TRayVector::iterator it;
    for( it = m_List.begin(); it != m_List.end(); ++it )
    {
      //TODO this should be a map with key defined as struct{name,position,type}
      // TODO this should be a map with key defined as ID, much better than struct{name,position,type} now that a ray has an ID
        // Can't add a line with a name + position + type that already exists in the list
        if( (*it).GetName() == r.GetName() && (*it).GetPosition() == r.GetPosition() && (*it).GetType() == r.GetType() )
	  throw GlobalException(INTERNAL_ERROR,Formatter()<<"Ray with name " << r.GetName() << " position " << r.GetPosition() << " and type " << r.GetType() << " already exists");
    }

    m_List.push_back( r );

}



void CRayCatalog::AddRayFromParams(const std::string& name,
				   const Float64& position,
				   const std::string& type,
				   const std::string& force,
				   const std::string& profileName,
				   const TAsymParams& asymParams,
				   const std::string& groupName,
				   const Float64& nominalAmplitude,
				   const std::string& velocityGroup,
				   const Float64& velocityOffset,
				   const bool& enableVelocityFit,
				   const Int32& id)
{
  int etype = -1;
    if (type == "E") etype=2;
    else if(type == "A") etype=1;
    else throw GlobalException(INTERNAL_ERROR, Formatter()<< "Bad ray type, should be in {A,E} : "<<type);

    int eforce = -1;
    if (force == "W") eforce = 1;
    else if (force == "S") eforce = 2;
    else throw GlobalException(INTERNAL_ERROR, Formatter()<< "Bad ray force, should be in {S,W} : "<<force);
    TAsymParams _asymParams = {1., 4.5, 0.};
    TAsymParams _asymFitParams = {2., 2., 0.};
    std::shared_ptr<CLineProfile> profile;

    if (profileName.find("ASYMFIXED") != std::string::npos) {
      profile = std::make_shared<CLineProfileASYM>(m_nSigmaSupport, asymParams, "mean");
    }
    else if (profileName == "SYM")
      profile = std::make_shared<CLineProfileSYM>(m_nSigmaSupport);
    else if (profileName == "LOR")
      profile = std::make_shared<CLineProfileLOR>(m_nSigmaSupport);
    else if (profileName == "ASYM"){

      profile = std::make_shared<CLineProfileASYM>(m_nSigmaSupport, _asymParams, "none");
    }
    else if (profileName == "ASYMFIT"){
      profile = std::make_shared<CLineProfileASYMFIT>(m_nSigmaSupport, _asymFitParams, "mean");
    }else{
      throw GlobalException(INTERNAL_ERROR, Formatter()<<"CRayCatalog::Load: Profile name "<<profileName<<" is no recognized.");
    }
    
      Add( CRay(name, position, etype, profile, eforce, velocityOffset, enableVelocityFit, groupName, nominalAmplitude, velocityGroup, id) );


  }


void CRayCatalog::Sort()
{
    sort(m_List.begin(), m_List.end());
}

void CRayCatalog::setLineAmplitude(const std::string& name,const Float64& nominalAmplitude)
{
    TRayVector::iterator it;
    for( it = m_List.begin(); it != m_List.end(); ++it )
    {
      if(it->GetName() == name) return it->setNominalAmplitude(nominalAmplitude);
    }
    throw GlobalException(INTERNAL_ERROR,Formatter()<<" Line " << name << " does not exist in catalog");
}
