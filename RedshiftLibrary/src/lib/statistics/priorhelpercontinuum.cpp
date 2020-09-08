#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/statistics/priorhelpercontinuum.h>

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

#include <RedshiftLibrary/log/log.h>

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
            Log.LogWarning("    CPriorHelperContinuum: rootFolder path does not exist: %s", rootFolder.string().c_str());
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
            Log.LogError("    CPriorHelperContinuum: AgaussMean path does not exist: %s", agaussfpath.string().c_str());
            throw std::runtime_error("    CPriorHelperContinuum: AgaussMean path does not exist");
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
            Log.LogError("    CPriorHelperContinuum: AgaussSigma path does not exist: %s", agaussfpath.string().c_str());
            throw std::runtime_error("    CPriorHelperContinuum: AgaussSigma path does not exist");
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
        Log.LogError("    CPriorHelperContinuum: SetTNameData failed for k=%d", k);
        throw std::runtime_error("    CPriorHelperContinuum: set bad tname index");
    }
    m_tplnames[k] = tname;
    return true;
}

bool CPriorHelperContinuum::SetEZTData(UInt32 k, std::vector<std::vector<Float64>> ezt_data)
{
    if(k>=m_data.size())
    {
        Log.LogError("    CPriorHelperContinuum: SetEZTData failed for k=%d", k);
        throw std::runtime_error("    CPriorHelperContinuum: set bad data index");
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
        Log.LogError("    CPriorHelperContinuum: SetAgaussmeanData failed for k=%d", k);
        throw std::runtime_error("    CPriorHelperContinuum: set bad data index");
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
        Log.LogError("    CPriorHelperContinuum: SetAgausssigmaData failed for k=%d", k);
        throw std::runtime_error("    CPriorHelperContinuum: set bad data index");
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
                    Log.LogError("    CPriorHelperContinuum: read n=%d cols, instead of %d", lineVals.size(), m_nEbv);
                    throw std::runtime_error("    CPriorHelperContinuum: read bad number of cols");
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



