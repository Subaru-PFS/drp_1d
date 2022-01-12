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


CRayCatalog::CRayCatalog()
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
                    filteredList.push_back(m_List[i].clone());
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


bool CRayCatalog::Add( const CRay& r )
{
    TRayVector::iterator it;
    for( it = m_List.begin(); it != m_List.end(); ++it )
    {
      //TODO this should be a map with key defined as struct{name,position,type}
      // TODO this should be a map with key defined as ID, much better than struct{name,position,type} now that a ray has an ID
        // Can't add a line with a name + position + type that already exists in the list
        if( (*it).GetName() == r.GetName() && (*it).GetPosition() == r.GetPosition() && (*it).GetType() == r.GetType() )
            return false;
    }

    m_List.push_back( r );

    return true;
}

/**
 * @brief parse_asymfixed
 * Parse a ASYMFIXED_XX_YY_ZZ profile string
 */
//static void parse_asymfixed(std::string &profileString, Float64 &width, Float64 &alpha, Float64 &delta)
static void parse_asymfixed(std::string &profileString, TAsymParams &asymParams)
{
    std::vector< Float64 > numbers;
    std::string temp;
    size_t start = 0, end;

    end = profileString.find("_", start);
    if (profileString.substr(start, end-start+1) != "ASYMFIXED_") {
        throw GlobalException(INTERNAL_ERROR,"XXX");
    }
    start = end + 1;
    end = profileString.find("_", start);
    temp = profileString.substr(start, end-start);
    asymParams.sigma = std::stod(temp);

    start = end + 1;
    end = profileString.find("_", start);
    temp = profileString.substr(start, end-start);
    asymParams.alpha = std::stod(temp);

    start = end + 1;
    end = profileString.length();
    temp = profileString.substr(start, end-start);
    asymParams.delta = std::stod(temp);
}

/**
 * @brief CRayCatalog::Load
 * Loads a line catalog in TSV format as of v0.4
 * @param filePath
 */
void CRayCatalog::Load( const char* filePath, Float64 nsigmasupport)
{
    std::ifstream file;

    // Clear current line list
    m_List.clear();

    if ( !exists( filePath ) ) {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Can't load line catalog : "<< filePath<<" does not exist.");
    }

    file.open( filePath, std::ifstream::in );
    if( file.rdstate() & ios_base::failbit ) {
      throw GlobalException(INTERNAL_ERROR,"file cannot be opened");
    }

    string line;

    Float64 ver = -1.0;
    // Read file line by line
    while( getline( file, line ) )
    {
        // manage version
        if(ver==-1.0){
            if(line.find("#version:0.4.0")!= std::string::npos)
            {
                ver = 0.4;
            }else if(line.find("#version:0.3.0")!= std::string::npos)
            {
                ver = 0.3;
            }else if(line.find("#version:")!= std::string::npos)
            {
                Log.LogDebug("Loading Line Catalog: unable to parse version line: %s", line.c_str());
            }
            continue;
        }

        if(ver!=0.4)
        {
	  throw GlobalException(INTERNAL_ERROR,Formatter()<<"Line catalog version (found ver="<<ver<<") is not supported");
        }

        // remove comments
        if(line.compare(0,1,"#",1)==0){
            continue;
        }
        char_separator<char> sep("\t");

        // Tokenize each line
        typedef tokenizer< char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );

        // Check if it's not a comment
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() && *it != "#" )
        {
            // Parse position
            double pos = 0.0;
            try{
                pos = lexical_cast<double>(*it);
            }catch (const bad_lexical_cast& e)
            {
	      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Bad file format : "<<filePath <<" : "<< e.what());
            }

            // Parse name
            ++it;
            string name;
            if( it != tok.end() )
                name = *it;
            else{
	      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Bad name in : "<< filePath);
            }

            // Parse type
            int Etype = 0;
            ++it;
            string type = "None";
            if( it != tok.end() )
                type = *it;
            if( strcmp(type.c_str(),"A")==0 )
                Etype = 1;
            else if( strcmp(type.c_str(),"E")==0 )
                Etype = 2;

            // Parse weak or strong
            int Eforce = 0;
            ++it;
            string strong = "None";
            if( it != tok.end() )
                strong = *it;
            if( strcmp(strong.c_str(),"W")==0 )
                Eforce = 1;
            else if( strcmp(strong.c_str(),"S")==0 )
                Eforce = 2;
            else{
	      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Bad force in :" << filePath);
            }

            std::string profileName = "SYM";
            TAsymParams asymParams = {NAN, NAN, NAN};
            std::string groupName = "-1";
            Float64 nominalAmplitude = 1.0;
            std::string velGroupName = "-1";

            //tmp: default values 
            TAsymParams _asymParams = {1., 4.5, 0.};
            TAsymParams _asymFitParams = {2., 2., 0.};


            std::shared_ptr<CLineProfile> profile;
            // Parse profile name
            ++it;
            if( it != tok.end() ){
                profileName = *it;
                if (profileName.find("ASYMFIXED") != std::string::npos) {
                    parse_asymfixed(profileName, asymParams);//reading params from catalog files
                    profile = std::make_shared<CLineProfileASYM>(nsigmasupport, asymParams, "mean");
                }
                else if (profileName == "SYM")
                    profile = std::make_shared<CLineProfileSYM>(nsigmasupport);
                else if (profileName == "LOR")
                    profile = std::make_shared<CLineProfileLOR>(nsigmasupport);
                else if (profileName == "ASYM"){
                    asymParams =  _asymParams;
                    profile = std::make_shared<CLineProfileASYM>(nsigmasupport, _asymParams, "none");
                }
                else if (profileName == "ASYMFIT"){
                    asymParams =  _asymFitParams; //using default values
                    profile = std::make_shared<CLineProfileASYMFIT>(nsigmasupport, _asymFitParams, "mean");
                }else{
                    throw GlobalException(INTERNAL_ERROR, Formatter()<<"CRayCatalog::Load: Profile name "<<profileName<<" is no recognized.");
                }
            }

            // Parse group name
            ++it;
            if( it != tok.end() )
            {
                groupName = *it;
                // Parse group line nominal amplitude
                ++it;
                if( it != tok.end() )
                {
                    try{
                        nominalAmplitude = lexical_cast<double>(*it);
                    }catch (bad_lexical_cast&){
                        Log.LogError( "Unable to read nominal amplitude value from file, setting as default (1.0)." );
                        nominalAmplitude = 1.0;
                    }
                }
                if(groupName=="" || groupName=="-1")
                    nominalAmplitude = 1.0;
            }

            // Parse velocity group name
            ++it;
            if( it != tok.end() )
                velGroupName = *it;

	    ++it;
	    int id=-1;
	    if( it != tok.end() ){
	        try{
		        id = lexical_cast<int>(*it);
            }catch (const bad_lexical_cast& e)
            {
                throw GlobalException(INTERNAL_ERROR,"Bad file format");
            }
        }
	    
        if( nominalAmplitude>0.0 ) //do not load a line with nominal amplitude = ZERO
	        Add( CRay(name, pos, Etype, profile, Eforce, -1, -1, -1, -1, -1, -1, groupName, nominalAmplitude, velGroupName, id) );
        }
    }
    file.close();

    if(ver<0.0)
    {
        throw GlobalException(INTERNAL_ERROR,"Invalid line catalog file (found ver");
    }
}

bool CRayCatalog::Save( const char* filePath )
{
    std::ofstream file;
    file.open( filePath, std::ofstream::out );
    if( file.rdstate() & ios_base::failbit )
    {
        return false;
    }

    file << "#version:0.4.0" << std::endl;
    file << "#lambda" << "\t" << "name" << "\t" << "type" << "\t" << "force" << "\t" << "profile" << "\t" << "amp_group" << "\t" << "nominal_ampl" << "\t" << "vel_group" << std::endl;
    for( int i = 0; i< m_List.size(); i++ )
    {
        file << m_List[i].GetPosition() << "\t";
        file << m_List[i].GetName() << "\t";
        if(m_List[i].GetType() == 1)
        {
            file << "A" << "\t";
        }else
        {
            file << "E" << "\t";
        }

        if(m_List[i].GetForce() == 1)
        {
            file << "W" << "\t";
        }else
        {
            file << "S" << "\t";
        }
        file << m_List[i].GetName() << "\t";

        file << m_List[i].GetGroupName() << "\t";
        file << m_List[i].GetNominalAmplitude() << "\t";

        file << m_List[i].GetVelGroupName();

        file << std::endl;

    }
    file.close();
    return true;
}

void CRayCatalog::Sort()
{
    sort(m_List.begin(), m_List.end());
}

