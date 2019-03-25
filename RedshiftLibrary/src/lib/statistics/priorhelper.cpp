#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/statistics/priorhelper.h>

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


CPriorHelper::CPriorHelper()
{
}

CPriorHelper::~CPriorHelper()
{

}

bool CPriorHelper::Init( std::string priorDirPath, Int32 type )
{
    m_type = type;

    bfs::path rootFolder( priorDirPath.c_str() );
    if(!bfs::exists(rootFolder))
    {
        Log.LogWarning("    CPriorHelper: rootFolder path does not exist: %s", rootFolder.string().c_str());
        Log.LogWarning("    CPriorHelper: priors won't be used");
        mInitFailed = true;
        return false;
    }

    std::vector<std::string> EZTfilesPathList;
    bfs::directory_iterator end_itr;
    std::string ezt_path = "";
    if(m_type==0)
    {
        ezt_path = "prior_continuum_hist_Ebmvc_Z";
    }else if(m_type==1)
    {
        ezt_path = "prior_lines_hist_Ebmvr_Z";
    }
    for ( bfs::directory_iterator itr( rootFolder/ezt_path.c_str() ); itr != end_itr; ++itr )
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

        std::string ezt_tag = "";
        std::string a_tag = "";
        std::string a_dirpath = "";
        if(m_type==0)
        {
            ezt_tag = "prior_continuum_hist_Ebmvc_Z_Tplc_";
            a_tag = "prior_continuum_gaussmean_Ac_Z_Tplc_";
            a_dirpath = "prior_continuum_gaussmean_Ac_Z";
        }else if(m_type==1)
        {
            ezt_tag = "prior_lines_hist_Ebmvr_Z_Tplr_";
            a_tag = "prior_lines_gaussmean_Ar_Z_Tplr_";
            a_dirpath = "prior_lines_gaussmean_Ar_Z";
        }

        boost::replace_all( fNameStr, ezt_tag.c_str(), a_tag.c_str());

        bfs::path agaussfpath = rootFolder/a_dirpath.c_str()/bfs::path(fNameStr);
        if(!bfs::exists(agaussfpath))
        {
            Log.LogError("    CPriorHelper: AgaussMean path does not exist: %s", agaussfpath.string().c_str());
            throw std::runtime_error("    CPriorHelper: AgaussMean path does not exist");
        }
        AGaussMeanfilesPathList.push_back(agaussfpath.string());
    }

    std::vector<std::string> AGaussSigmafilesPathList;
    for(UInt32 k=0; k<EZTfilesPathList.size(); k++)
    {
        bfs::path fPath = EZTfilesPathList[k];
        std::string fNameStr = fPath.filename().c_str();

        std::string ezt_tag = "";
        std::string a_tag = "";
        std::string a_dirpath = "";
        if(m_type==0)
        {
            ezt_tag = "prior_continuum_hist_Ebmvc_Z_Tplc_";
            a_tag = "prior_continuum_gausssigma_Ac_Z_Tplc_";
            a_dirpath = "prior_continuum_gausssigma_Ac_Z";
        }else if(m_type==1)
        {
            ezt_tag = "prior_lines_hist_Ebmvr_Z_Tplr_";
            a_tag = "prior_lines_gausssigma_Ar_Z_Tplr_";
            a_dirpath = "prior_lines_gausssigma_Ar_Z";
        }

        boost::replace_all( fNameStr, ezt_tag.c_str(), a_tag.c_str());

        bfs::path agaussfpath = rootFolder/a_dirpath.c_str()/bfs::path(fNameStr);
        if(!bfs::exists(agaussfpath))
        {
            Log.LogError("    CPriorHelper: AgaussSigma path does not exist: %s", agaussfpath.string().c_str());
            throw std::runtime_error("    CPriorHelper: AgaussSigma path does not exist");
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


    std::string z_dirpath = "";
    std::string z_filepath = "";
    if(m_type==0)
    {
        z_dirpath = "prior_continuum_hist_Z";
        z_filepath = "prior_continuum_hist_Z.txt";
    }else if(m_type==1)
    {
        z_dirpath = "prior_lines_hist_Z";
        z_filepath = "prior_lines_hist_Z.txt";
    }
    bfs::path pz_fpath = rootFolder/z_dirpath.c_str()/z_filepath.c_str();
    if(!bfs::exists(pz_fpath))
    {
        Log.LogError("    CPriorHelper: Pz path does not exist: %s", pz_fpath.string().c_str());
        throw std::runtime_error("    CPriorHelper: Pz path does not exist");
    }else{
        std::vector<Float64> read_buffer;
        bool ret = LoadFileZ(pz_fpath.string().c_str(), read_buffer);
        if(!ret)
        {
            Log.LogError("Unable to load the Pz prior data. aborting...");
            mInitFailed = true;
            return false;
        }
        SetPzData(read_buffer);
        mInitFailed = false;
    }

    return true;
}


bool CPriorHelper::SetBetaA(Float64 beta)
{
    m_betaA=beta;
    return true;
}

bool CPriorHelper::SetBetaTE(Float64 beta)
{
    m_betaTE=beta;
    return true;
}

bool CPriorHelper::SetBetaZ(Float64 beta)
{
    m_betaZ=beta;
    return true;
}

bool CPriorHelper::SetSize(UInt32 size)
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

    m_data_pz.clear();
    for(UInt32 kz=0; kz<m_nZ; kz++)
    {
        m_data_pz.push_back(0.0);
    }

    m_tplnames.clear();
    m_tplnames.resize(size);
    return true;
}


bool CPriorHelper::SetTNameData(UInt32 k, std::string tname)
{
    if(k>=m_tplnames.size())
    {
        Log.LogError("    CPriorHelper: SetTNameData failed for k=%d", k);
        throw std::runtime_error("    CPriorHelper: set bad tname index");
    }
    m_tplnames[k] = tname;
    return true;
}

bool CPriorHelper::SetEZTData(UInt32 k, std::vector<std::vector<Float64>> ezt_data)
{
    if(k>=m_data.size())
    {
        Log.LogError("    CPriorHelper: SetEZTData failed for k=%d", k);
        throw std::runtime_error("    CPriorHelper: set bad data index");
    }

    for(UInt32 kz=0; kz<m_nZ; kz++)
    {
        for(UInt32 ke=0; ke<m_nEbv; ke++)
        {
            m_data[k][kz][ke].priorTZE = ezt_data[kz][ke];
        }
    }

    return true;
}

bool CPriorHelper::SetAGaussmeanData(UInt32 k, std::vector<std::vector<Float64>> agaussmean_data)
{
    if(k>=m_data.size())
    {
        Log.LogError("    CPriorHelper: SetAgaussmeanData failed for k=%d", k);
        throw std::runtime_error("    CPriorHelper: set bad data index");
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

bool CPriorHelper::SetAGausssigmaData(UInt32 k, std::vector<std::vector<Float64>> agausssigma_data)
{
    if(k>=m_data.size())
    {
        Log.LogError("    CPriorHelper: SetAgausssigmaData failed for k=%d", k);
        throw std::runtime_error("    CPriorHelper: set bad data index");
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

bool CPriorHelper::SetPzData(std::vector<Float64> z_data)
{
    if(z_data.size()!=m_data_pz.size())
    {
        Log.LogError("    CPriorHelper: SetPzData failed for bad data size" );
        throw std::runtime_error("    CPriorHelper: set pz bad data size");
    }

    for(UInt32 kz=0; kz<m_nZ; kz++)
    {
        m_data_pz[kz] = z_data[kz];
    }

    return true;
}


bool CPriorHelper::LoadFileEZ( const char* filePath, std::vector<std::vector<Float64>>& data)
{
    bool verboseRead=false;
    Log.LogDetail("    CPriorHelper: start load prior file: %s", filePath);
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
                        Log.LogInfo("    CPriorHelper: read line=%d, col=%d : valf=%e", nlinesRead, lineVals.size()-1, x);
                    }
                }
                if(lineVals.size()!=m_nEbv)
                {
                    Log.LogError("    CPriorHelper: read n=%d cols, instead of %d", lineVals.size(), m_nEbv);
                    throw std::runtime_error("    CPriorHelper: read bad number of cols");
                }
                nlinesRead++;
                data.push_back(lineVals);
                if(verboseRead)
                {
                    Log.LogInfo("    CPriorHelper: read n=%d cols", lineVals.size());
                }
            }
        }
        file.close();
    }

    return loadSuccess;
}

bool CPriorHelper::LoadFileZ(const char* filePath , std::vector<Float64>& data)
{
    bool verboseRead=false;
    Log.LogDetail("    CPriorHelper: start load prior file: %s", filePath);
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

                std::istringstream iss( line );
                Float64 x;
                iss >> x;
                if(std::isnan(x) || std::isinf(x) || x!=x || x<0.)
                {
                    x = m_priorminval;
                }

                if(verboseRead)
                {
                    Log.LogInfo("    CPriorHelper: read line=%d : valf=%e", nlinesRead, x);
                }

                nlinesRead++;
                data.push_back(x);
            }
        }
        file.close();

        if(nlinesRead!=m_nZ)
        {
            Log.LogError("    CPriorHelper: read n=%d lines", nlinesRead);
            throw std::runtime_error("    CPriorHelper: read bad number of lines");
        }
    }
    return loadSuccess;
}


/**
 * @brief CPriorHelper::GetTplPriorData
 * @param tplname
 * @param redshifts
 * @param zePriorData
 * @param outsideZRangeExtensionMode: 0=extend 0 value for z<m_z0, and n-1 value for z>m_z0+m_dZ*m_nZ, 1=return error if z outside prior range
 * @return
 */
bool CPriorHelper::GetTplPriorData(std::string tplname,
                                            std::vector<Float64> redshifts,
                                            TPriorZEList& zePriorData,
                                            Int32 outsideZRangeExtensionMode)
{
    bool verbose=false;
    if(m_betaA<=0.0 && m_betaTE<=0.0 && m_betaZ<=0.0)
    {
        Log.LogError("    CPriorHelper: beta coeff all zero (betaA=%e, betaTE=%e, betaZ=%e)", m_betaA, m_betaTE, m_betaZ);
        zePriorData.clear();
        return false;
    }

    if(mInitFailed)
    {
        Log.LogDetail("    CPriorHelper: init. failed, unable to provide priors. Priors won't be used.");
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
        Log.LogError("    CPriorHelper: unable to match this tplname in priors names list : %s", tplname.c_str());
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
               Log.LogError("    CPriorHelper: unable to match this redshift in prior list : %e", redshifts[kz]);
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
           Log.LogInfo("    CPriorHelper: get prior for z=%f: found idz=%d", redshifts[kz], idz);
       }
       TPriorEList dataz = m_data[idx][idz];
       for(UInt32 icol=0; icol<m_nEbv; icol++)
       {
           if(verbose)
           {
               Log.LogInfo("    CPriorHelper: get prior for tpl=%s", tplname.c_str());
               Log.LogInfo("    CPriorHelper: get prior idTpl=%d, idz=%d, idebmv=%d : valf=%e", idx, idz, icol, dataz[icol].priorTZE);
           }

           Float64 logPTE = 0.;
           Float64 logPA = 0.;
           Float64 logPZ = 0.;
           if(dataz[icol].A_sigma>0.0)
           {
               logPA =-0.5*log(2*M_PI) - log(dataz[icol].A_sigma);
           }
           if(dataz[icol].priorTZE>0.0)
           {
               logPTE = log(dataz[icol].priorTZE);

               if(std::isnan(logPTE) || logPTE!=logPTE || std::isinf(logPTE))
               {
                   Log.LogError("    CPriorHelper: logP_TZE is NAN (priorTZE=%e, logP_TZE=%e)",
                                dataz[icol].priorTZE,
                                logPTE);
                   throw std::runtime_error("    CPriorHelper: logP_TZE is NAN or inf, or invalid");
               }
           }

           if(m_data_pz[idz]>0.0)
           {
               logPZ = log(m_data_pz[idz]/m_dz);
           }

           dataz[icol].logprior_precompA = logPA;
           dataz[icol].logprior_precompTE = logPTE;
           dataz[icol].logprior_precompZ = logPZ;
           dataz[icol].betaA = m_betaA;
           dataz[icol].betaTE = m_betaTE;
           dataz[icol].betaZ = m_betaZ;
       }
       zePriorData.push_back(dataz);
    }


    return true;
}


bool CPriorHelper::GetTZEPriorData(std::string tplname,
                                   Int32 EBVIndexfilter,
                                   Float64 redshift,
                                   SPriorTZE& tzePrioData,
                                   Int32 outsideZRangeExtensionMode)
{
    if(EBVIndexfilter<0 || EBVIndexfilter>m_nEbv-1)
    {
        Log.LogError("    CPriorHelper: Bad EBV index requested =%d (nEBV=%d)", EBVIndexfilter, m_nEbv);
        throw std::runtime_error("    CPriorHelper: Bad EBV index requested");
    }
    std::vector<Float64> redshifts(1, redshift);
    TPriorZEList zePriorData;
    GetTplPriorData(tplname, redshifts, zePriorData, outsideZRangeExtensionMode);

    tzePrioData = zePriorData[0][EBVIndexfilter];

    return true;
}


