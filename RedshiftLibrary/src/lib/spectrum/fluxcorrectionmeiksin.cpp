#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"

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
#include <numeric>
namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;
using namespace boost;


CSpectrumFluxCorrectionMeiksin::CSpectrumFluxCorrectionMeiksin()
{
    m_LambdaMin = 200.0;
    m_LambdaMax = 1299.0;
}

CSpectrumFluxCorrectionMeiksin::~CSpectrumFluxCorrectionMeiksin()
{

}

Bool CSpectrumFluxCorrectionMeiksin::Init( std::string calibrationPath, const std::shared_ptr<const CLSF>& lsf )
{

    bfs::path calibrationFolder( calibrationPath.c_str() );
    std::vector<std::string> fileNamesList;
    fileNamesList.push_back("Meiksin_Var_curves_2.0.txt");
    fileNamesList.push_back("Meiksin_Var_curves_2.5.txt");
    fileNamesList.push_back("Meiksin_Var_curves_3.0.txt");
    fileNamesList.push_back("Meiksin_Var_curves_3.5.txt");
    fileNamesList.push_back("Meiksin_Var_curves_4.0.txt");
    fileNamesList.push_back("Meiksin_Var_curves_4.5.txt");
    fileNamesList.push_back("Meiksin_Var_curves_5.0.txt");
    fileNamesList.push_back("Meiksin_Var_curves_5.5.txt");
    fileNamesList.push_back("Meiksin_Var_curves_6.0.txt");
    fileNamesList.push_back("Meiksin_Var_curves_6.5.txt");
    fileNamesList.push_back("Meiksin_Var_curves_7.0.txt");

    for(UInt32 k=0; k<fileNamesList.size(); k++)
    {
        bfs::path fPath = ( calibrationFolder/"igm"/"IGM_variation_curves_meiksin"/fileNamesList[k].c_str() ).string();
        std::string fPathStr = (fPath).string();

        bool ret = LoadCurvesinIncreasingExtinctionOrder(fPathStr.c_str());
        if(!ret)
        {
            Log.LogError("Unable to load the Meiksin flux correction data. aborting...");
            meiksinInitFailed = true;
            return false;
        }
        meiksinInitFailed = false;
    }

    ConvolveAll(lsf);
    return true;
}

/**
 * Important: igm curves should be loaded in the increasing ordre of their extinction per bin of z,
 * i.e., from the least extinction curve to the highest extinction curve 
 *         
 * */
Bool CSpectrumFluxCorrectionMeiksin::LoadCurvesinIncreasingExtinctionOrder( const char* filePath )
{
    bool loadSuccess=true;
    std::ifstream file;
    file.open( filePath, std::ifstream::in );
    bool fileOpenFailed = file.rdstate() & std::ios_base::failbit;
    if(fileOpenFailed)
    {
        loadSuccess = false;
    }else
    {
        MeiksinCorrection _correction;
        TFloat64List fluxcorr0;
        TFloat64List fluxcorr1;
        TFloat64List fluxcorr2;
        TFloat64List fluxcorr3;
        TFloat64List fluxcorr4;
        TFloat64List fluxcorr5;
        TFloat64List fluxcorr6;

        std::string line;
        // Read file line by line
        while( getline( file, line ) )
        {
            if( !boost::starts_with( line, "#" ) )
            {
                std::istringstream iss( line );
                Float64 x, y0, y1, y2, y3, y4, y5, y6;
                iss >> x >> y0 >> y1 >> y2 >> y3 >> y4 >> y5 >> y6;
                _correction.lbda.push_back(x);
                fluxcorr0.push_back(y0);
                fluxcorr1.push_back(y1);
                fluxcorr2.push_back(y2);
                fluxcorr3.push_back(y3);
                fluxcorr4.push_back(y4);
                fluxcorr5.push_back(y5);
                fluxcorr6.push_back(y6);
            }
        }
        file.close();
        _correction.fluxcorr.push_back(fluxcorr0);
        _correction.fluxcorr.push_back(fluxcorr1);
        _correction.fluxcorr.push_back(fluxcorr2);
        _correction.fluxcorr.push_back(fluxcorr3);
        _correction.fluxcorr.push_back(fluxcorr4);
        _correction.fluxcorr.push_back(fluxcorr5);
        _correction.fluxcorr.push_back(fluxcorr6);

        m_rawCorrections.push_back(_correction);
    }

    return loadSuccess;
}

/**
 * @brief CSpectrumFluxCorrectionMeiksin::GetRedshiftIndex
 * @param z
 *
 * Returns the index corresponding to the table to be used for that redshift value
 *
 * @return
 */
Int32 CSpectrumFluxCorrectionMeiksin::GetRedshiftIndex(Float64 z) const
{
    Int32 index = -1;

    Float64 zStart = 2.0;
    Float64 zStep = 0.5;
    Float64 zStop = 7.0;
    if(z<zStart){
        index = 0;
    }else if( z>=zStart && z<zStop)
    {
        index = (Int32)( (z-zStart)/zStep + 1.0);
    }else if( z>=zStop )
    {
        index = (Int32)((zStop-zStart)/zStep);
    }
    return index;
}

//TBC once validated
TFloat64List CSpectrumFluxCorrectionMeiksin::GetLSFProfileVector(Float64 lambda0_rest, Float64 z_bin_meiksin, const std::shared_ptr<const CLSF>& lsf )
{
    if(!lsf->IsValid()){
        Log.LogError("LSF is not valid");
        throw std::runtime_error("LSF is not valid");
    }
    Float64 lambda0_obs = lambda0_rest*(1+z_bin_meiksin);
    Float64 sigma_obs = lsf->GetWidth(lambda0_obs);
    Float64 sigmaSupport = lsf->GetProfile()->GetNSigmaSupport();

    Float64 lbdastep_rest = 1.;//value in angstrom based on calibration-igm files
    Int32 Nhalf = std::round(sigmaSupport*sigma_obs/(1+z_bin_meiksin)/lbdastep_rest);
    Int32 len = 2*Nhalf+1;

    //change to observedframe
    TFloat64List lambdas_obs(len);
    for(Int32 i=0; i<len; i++)
    {
        lambdas_obs[i] = (lambda0_rest +(i-Nhalf)*lbdastep_rest)*(1+z_bin_meiksin);
    }
    Float64 norm = 0, v; 
    m_kernel.resize(len);
    for(Int32 i = 0; i<len; i++)
    {
        //getLineProfile expects lbda in observedframe
        v = lsf->GetLineProfile(lambdas_obs[i], lambda0_obs, sigma_obs);
        m_kernel[i] = v;
        norm+= v; 
    }
    //normalizing  values
    Float64 inv_norm = 1/norm;
    for(Int32 i =0; i<len; i++){
        m_kernel[i] *= inv_norm;
    }
    return m_kernel;
}
/**
 * f*g[0] = f[0].g[0]
 * f*g[1] = f[0].g[1] + f[1].g[0] = sum(f[i].g[j]) for i, j belonging the intersection
*/
TFloat64List CSpectrumFluxCorrectionMeiksin::Convolve(const TFloat64List& arr, const TFloat64List& kernel)
{
    if(!arr.size() || !kernel.size())
    {
        throw std::runtime_error("Cannot convolve: either kernel or array is empty. ");
    }
    Int32 n = arr.size(), Nhalf = int(kernel.size()/2);
    TFloat64List convolvedArr(n);
    
    Float64 tmp;
    for(Int32 i = 0; i < n; i++){
        tmp = 0.0;//sum over the intersection area
        for(Int32 j = -Nhalf; j<=Nhalf; j++)//center kernel at arr[i]
        {
            if(i+j>=0)
            {
                if(i+j>=n){
                    tmp += kernel[Nhalf+j];
                }else{
                    tmp += arr[i+j]*kernel[Nhalf+j];
                }
            }
        }
        convolvedArr[i] = tmp;
    }
    return convolvedArr;
}
//loop over lambda values
//for each lambda0, compute kernel_lambda0 and then multiply it by the igm curve
TFloat64List CSpectrumFluxCorrectionMeiksin::ApplyAdaptativeKernel(const TFloat64List& arr, 
                                                                 const Float64 z_center, 
                                                                 const std::shared_ptr<const CLSF>& lsf,
                                                                 const TFloat64List& lambdas)
{
    if(!arr.size()){
        Log.LogError("Cannot convolve: either kernel or array is empty. ");
        throw std::runtime_error("Cannot convolve: either kernel or array is empty. ");
    }

    Int32 n = arr.size(), Nhalf=-1;
    TFloat64List convolvedArr(n);
    
    Float64 tmp;
    for(Int32 i = 0; i < n; i++){
        Float64 lambda0 = lambdas[i];//lambda restframe
        //compute the adpative kernel at lambda0
        GetLSFProfileVector(lambda0, z_center, lsf);//resulting kernel saved in m_kernel
        Nhalf = int(m_kernel.size()/2);
        if(!m_kernel.size()){
            Log.LogError("Cannot convolve: either kernel or array is empty. ");  
            throw std::runtime_error("Cannot convolve: either kernel or array is empty. ");
        }

        tmp = 0.0;//sum over the intersection area
        for(Int32 j = -Nhalf; j<=Nhalf; j++)//center kernel at arr[j]
        {
            if(i+j>=0)
            { 
                if(i+j>=n){
                    tmp += m_kernel[Nhalf+j];
                }else{
                    tmp += arr[i+j]*m_kernel[Nhalf+j];
                }
            }else{
                tmp += arr[0]*m_kernel[Nhalf+j];
            }
        }
        convolvedArr[i] = tmp;
    }
    return convolvedArr;
}
void CSpectrumFluxCorrectionMeiksin::ConvolveAll(const std::shared_ptr<const CLSF>& lsf)
{
    //to be decided if we create two objects or we keep one m_correction
    m_corrections.resize(m_rawCorrections.size());

    TFloat64List meiksin_Bins = GetSegmentsStartRedshiftList();
    Float64 zstep = meiksin_Bins[2]-meiksin_Bins[1];

    //iterate over the redshift list
    Float64 z_center;
    for(Int32 i = 0; i<m_rawCorrections.size(); i++)
    {
        if(i==0) z_center = meiksin_Bins[i+1]- zstep/2.;
        else     z_center = meiksin_Bins[i]  + zstep/2.;
        
        m_corrections[i].lbda = m_rawCorrections[i].lbda;
        
        for(Int32 j = 0; j<m_rawCorrections[i].fluxcorr.size(); j++) //iterating over the different curves
        {   
            TFloat64List a = ApplyAdaptativeKernel(m_rawCorrections[i].fluxcorr[j], z_center, lsf, m_corrections[i].lbda);
            m_corrections[i].fluxcorr.push_back(std::move(a));
        }
    }
    return;
}