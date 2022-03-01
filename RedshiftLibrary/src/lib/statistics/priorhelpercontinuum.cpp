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
#include "RedshiftLibrary/statistics/priorhelpercontinuum.h"
#include "RedshiftLibrary/common/flag.h"

#include <algorithm>    // std::sort
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "RedshiftLibrary/log/log.h"

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;
using namespace boost;


CPriorHelperContinuum::CPriorHelperContinuum()
{
}

CPriorHelperContinuum::~CPriorHelperContinuum()
{

}

bool CPriorHelperContinuum::Init( std::string priorDirPath )
{
    bfs::path rootFolder( priorDirPath.c_str() );
    if(!bfs::exists(rootFolder))
    {
        if(!rootFolder.string().empty())
        {
            Flag.warning(Flag.INVALID_FOLDER_PATH, Formatter()<<"    CPriorHelperContinuum::"<<__func__<<": rootFolder path does not exist: " << rootFolder.string().c_str());
            Log.LogDetail("    CPriorHelperContinuum: priors won't be used");
        }
        mInitFailed = true;
        return false;
    }

    std::vector<std::string> EZTfilesPathList;
    bfs::directory_iterator end_itr;
    for ( bfs::directory_iterator itr( rootFolder/"prior_continuum_hist_Ebmvc_Z" ); itr != end_itr; ++itr )
    {
        if ( !is_directory( itr->status() ) )
        {
            EZTfilesPathList.push_back(itr->path().c_str());
        }
    }

    std::vector<std::string> AGaussMeanfilesPathList;
    for(UInt32 k=0; k<EZTfilesPathList.size(); k++)
    {
        bfs::path fPath = EZTfilesPathList[k];
        std::string fNameStr = fPath.filename().c_str();
        boost::replace_all( fNameStr, "prior_continuum_hist_Ebmvc_Z_Tplc_", "prior_continuum_gaussmean_Ac_Z_Tplc_");

        bfs::path agaussfpath = rootFolder/"prior_continuum_gaussmean_Ac_Z"/bfs::path(fNameStr);
        if(!bfs::exists(agaussfpath))
        {
	  throw GlobalException(INTERNAL_ERROR,Formatter()<<"CPriorHelperContinuum: AgaussMean path does not exist: "<< agaussfpath.string());
        }
        AGaussMeanfilesPathList.push_back(agaussfpath.string());
    }

    std::vector<std::string> AGaussSigmafilesPathList;
    for(UInt32 k=0; k<EZTfilesPathList.size(); k++)
    {
        bfs::path fPath = EZTfilesPathList[k];
        std::string fNameStr = fPath.filename().c_str();
        boost::replace_all( fNameStr, "prior_continuum_hist_Ebmvc_Z_Tplc_", "prior_continuum_gausssigma_Ac_Z_Tplc_");

        bfs::path agaussfpath = rootFolder/"prior_continuum_gausssigma_Ac_Z"/bfs::path(fNameStr);
        if(!bfs::exists(agaussfpath))
        {
	  throw GlobalException(INTERNAL_ERROR,Formatter()<<"CPriorHelperContinuum: AgaussSigma path does not exist: "<<agaussfpath.string());
        }
        AGaussSigmafilesPathList.push_back(agaussfpath.string());
    }

    //allocate the data buffer
    SetSize(EZTfilesPathList.size());


    //set the template names
    for(UInt32 k=0; k<EZTfilesPathList.size(); k++)
    {
        bfs::path fPath = EZTfilesPathList[k];
        std::string fNameStr = fPath.filename().c_str();
        SetTNameData(k, fNameStr);
    }

    //read the EZT data from files
    for(UInt32 k=0; k<EZTfilesPathList.size(); k++)
    {
        bfs::path fPath = EZTfilesPathList[k];
        std::string fPathStr = (fPath).string();

        std::vector<std::vector<Float64>> read_buffer;
        bool ret = LoadFileEZ(fPathStr.c_str(), read_buffer);
        if(!ret)
        {
            Log.LogError("Unable to load the EZT continuum prior data. aborting...");
            mInitFailed = true;
            return false;
        }
        SetEZTData(k, read_buffer);
        mInitFailed = false;
    }


    //read the AGaussMean data from files
    for(UInt32 k=0; k<AGaussMeanfilesPathList.size(); k++)
    {
        bfs::path fPath = AGaussMeanfilesPathList[k];
        std::string fPathStr = (fPath).string();

        std::vector<std::vector<Float64>> read_buffer;
        bool ret = LoadFileEZ(fPathStr.c_str(), read_buffer);
        if(!ret)
        {
            Log.LogError("Unable to load the A-gaussmean continuum prior data. aborting...");
            mInitFailed = true;
            return false;
        }
        SetAGaussmeanData(k, read_buffer);
        mInitFailed = false;
    }

    //read the AGaussSigma data from files
    for(UInt32 k=0; k<AGaussSigmafilesPathList.size(); k++)
    {
        bfs::path fPath = AGaussSigmafilesPathList[k];
        std::string fPathStr = (fPath).string();

        std::vector<std::vector<Float64>> read_buffer;
        bool ret = LoadFileEZ(fPathStr.c_str(), read_buffer);
        if(!ret)
        {
            Log.LogError("Unable to load the A-gausssigma continuum prior data. aborting...");
            mInitFailed = true;
            return false;
        }
        SetAGausssigmaData(k, read_buffer);
        mInitFailed = false;
    }
    return true;
}


bool CPriorHelperContinuum::SetBeta(Float64 beta)
{
    m_beta=beta;
    return true;
}

bool CPriorHelperContinuum::SetSize(UInt32 size)
{
    m_data.clear();
    for(UInt32 k=0; k<size; k++)
    {
        TPriorZEList _zelist;
        for(UInt32 kz=0; kz<m_nZ; kz++)
        {
            std::vector<SPriorTZE> _elist;
            for(UInt32 ke=0; ke<m_nEbv; ke++)
            {
                SPriorTZE _tze;
                _elist.push_back(_tze);
            }
            _zelist.push_back(_elist);
        }
        m_data.push_back(_zelist);
    }

    m_tplnames.clear();
    m_tplnames.resize(size);
    return true;
}


bool CPriorHelperContinuum::SetTNameData(UInt32 k, std::string tname)
{
    if(k>=m_tplnames.size())
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"CPriorHelperContinuum: SetTNameData failed for k="<< k);
    }
    m_tplnames[k] = tname;
    return true;
}

bool CPriorHelperContinuum::SetEZTData(UInt32 k, std::vector<std::vector<Float64>> ezt_data)
{
    if(k>=m_data.size())
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"CPriorHelperContinuum: SetEZTData failed for k="<<k);
    }

    for(UInt32 kz=0; kz<m_nZ; kz++)
    {
        for(UInt32 ke=0; ke<m_nEbv; ke++)
        {
            m_data[k][kz][ke].logpriorTZE = ezt_data[kz][ke];
        }
    }

    return true;
}

bool CPriorHelperContinuum::SetAGaussmeanData(UInt32 k, std::vector<std::vector<Float64>> agaussmean_data)
{
    if(k>=m_data.size())
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"CPriorHelperContinuum: SetAgaussmeanData failed for k="<< k);
    }

    for(UInt32 kz=0; kz<m_nZ; kz++)
    {
        for(UInt32 ke=0; ke<m_nEbv; ke++)
        {
            m_data[k][kz][ke].A_mean = agaussmean_data[kz][ke];
        }
    }

    return true;
}

bool CPriorHelperContinuum::SetAGausssigmaData(UInt32 k, std::vector<std::vector<Float64>> agausssigma_data)
{
    if(k>=m_data.size())
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"CPriorHelperContinuum: SetAgausssigmaData failed for k="<<k);
    }

    for(UInt32 kz=0; kz<m_nZ; kz++)
    {
        for(UInt32 ke=0; ke<m_nEbv; ke++)
        {
            m_data[k][kz][ke].A_sigma = agausssigma_data[kz][ke];
        }
    }

    return true;
}

bool CPriorHelperContinuum::LoadFileEZ( const char* filePath, std::vector<std::vector<Float64>>& data)
{
    bool verboseRead=false;
    Log.LogDetail("    CPriorHelperContinuum: start load prior file: %s", filePath);
    bool loadSuccess=true;
    std::ifstream file;
    file.open( filePath, std::ifstream::in );
    bool fileOpenFailed = file.rdstate() & std::ios_base::failbit;
    if(fileOpenFailed)
    {
        loadSuccess = false;
    }else
    {
        Int32 nlinesRead = 0;
        std::string line;
        // Read file line by line
        while( getline( file, line ) )
        {
            if( !boost::starts_with( line, "#" ) )
            {

                std::vector<Float64> lineVals;
                std::istringstream iss( line );
                for(UInt32 icol=0; icol<m_nEbv; icol++)
                {
                    Float64 x;
                    iss >> x;
                    if(std::isnan(x) || std::isinf(x) || x!=x || x<=0.)
                    {
                        x = m_priorminval;
                    }
                    lineVals.push_back(x);

                    if(verboseRead)
                    {
                        Log.LogDetail("    CPriorHelperContinuum: read line=%d, col=%d : valf=%e", nlinesRead, lineVals.size()-1, x);
                    }
                }
                if(lineVals.size()!=m_nEbv)
                {
		  throw GlobalException(INTERNAL_ERROR,Formatter()<<"CPriorHelperContinuum: read n="<< lineVals.size()<<" cols, instead of "<< m_nEbv);
                }
                nlinesRead++;
                data.push_back(lineVals);
                if(verboseRead)
                {
                    Log.LogDetail("    CPriorHelperContinuum: read n=%d cols", lineVals.size());
                }
            }
        }
        file.close();
    }

    return loadSuccess;
}

/**
 * @brief CPriorHelperContinuum::GetTplPriorData
 * @param tplname
 * @param redshifts
 * @param zePriorData
 * @param outsideZRangeExtensionMode: 0=extend 0 value for z<m_z0, and n-1 value for z>m_z0+m_dZ*m_nZ, 1=return error if z outside prior range
 * @return
 */
bool CPriorHelperContinuum::GetTplPriorData(std::string tplname,
                                            std::vector<Float64> redshifts,
                                            TPriorZEList& zePriorData,
                                            Int32 outsideZRangeExtensionMode)
{
    bool verbose=false;
    if(m_beta<=0.0)
    {
        Log.LogError("    CPriorHelperContinuum: beta coeff not initialized correctly (=%e)", m_beta);
        zePriorData.clear();
        return false;
    }

    if(mInitFailed)
    {
        Log.LogDetail("    CPriorHelperContinuum: init. failed, unable to provide priors. Priors won't be used.");
        zePriorData.clear();
        return true;
    }

    //find idx for tplname
    UInt32 idx=-1;
    for(Int32 k=0; k<m_tplnames.size(); k++)
    {
        std::size_t foundstra = m_tplnames[k].find(tplname.c_str());
        if (foundstra!=std::string::npos){
            idx=k;
            break;
        }
    }
    if(idx<0 || idx>=m_tplnames.size())
    {
        Log.LogError("    CPriorHelperContinuum: unable to match this tplname in priors names list : %s", tplname.c_str());
        return false;
    }

    zePriorData.clear();
    Int32 idz=-1;
    for(Int32 kz=0; kz<redshifts.size(); kz++)
    {

       idz = Int32( (redshifts[kz]-m_z0)/m_dz );
       if(outsideZRangeExtensionMode==1)
       {
           if(idz<0 || idz>=m_nZ)
           {
               Log.LogError("    CPriorHelperContinuum: unable to match this redshift in prior list : %e", redshifts[kz]);
               return false;
           }
       }
       if(outsideZRangeExtensionMode==0)
       {
           if(idz<0)
           {
               idz=0;
           }
           if(idz>=m_nZ)
           {
               idz=m_nZ;
           }

       }
       if(verbose)
       {
           Log.LogDetail("    CPriorHelperContinuum: get prior for z=%f: found idz=%d", redshifts[kz], idz);
       }
       TPriorEList dataz = m_data[idx][idz];
       for(UInt32 icol=0; icol<m_nEbv; icol++)
       {
           if(verbose)
           {
               Log.LogDetail("    CPriorHelperContinuum: get prior for tpl=%s", tplname.c_str());
               Log.LogDetail("    CPriorHelperContinuum: get prior idTpl=%d, idz=%d, idebmv=%d : valf=%e", idx, idz, icol, dataz[icol].logpriorTZE);
           }
            dataz[icol].logpriorTZE /= m_dz;
            dataz[icol].logpriorTZE = m_beta*log(dataz[icol].logpriorTZE);
       }
       zePriorData.push_back(dataz);
    }


    return true;
}



