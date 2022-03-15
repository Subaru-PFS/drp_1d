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
#include "RedshiftLibrary/common/indexing.h"
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


CSpectrumFluxCorrectionMeiksin::CSpectrumFluxCorrectionMeiksin()//should be removed
{
    m_LambdaMin = 200.0;
    m_LambdaMax = 1299.0;
}

CSpectrumFluxCorrectionMeiksin::CSpectrumFluxCorrectionMeiksin(std::vector<MeiksinCorrection> meiksinCorrectionCorves):
m_rawCorrections(std::move(meiksinCorrectionCorves)),
m_LambdaMin(200.0),
m_LambdaMax(1299.0)
{}

void CSpectrumFluxCorrectionMeiksin::init(const std::shared_ptr<const CLSF>& lsf, TFloat64Range& convolRange)
{
    m_convolRange = convolRange;
    convolveAll(lsf);
}
/**
 * @brief CSpectrumFluxCorrectionMeiksin::getRedshiftIndex
 * @param z
 *
 * Returns the index corresponding to the table to be used for that redshift value
 *
 * @return
 */
Int32 CSpectrumFluxCorrectionMeiksin::getRedshiftIndex(Float64 z) const
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
TFloat64List CSpectrumFluxCorrectionMeiksin::getLSFProfileVector(Float64 lambda0_rest, Float64 z_bin_meiksin, const std::shared_ptr<const CLSF>& lsf )
{
    if(!lsf->IsValid()){
        throw GlobalException(INTERNAL_ERROR,"LSF is not valid");
    }
    Float64 lambda0_obs = lambda0_rest*(1+z_bin_meiksin);
    Float64 sigma_obs = lsf->GetWidth(lambda0_obs);
    Float64 sigmaSupport = lsf->GetProfile().GetNSigmaSupport();

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
TFloat64List CSpectrumFluxCorrectionMeiksin::convolve(const TFloat64List& arr, const TFloat64List& kernel)
{
    if(!arr.size() || !kernel.size())
    {
        throw GlobalException(INTERNAL_ERROR,"Cannot convolve: either kernel or array is empty. ");
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
TFloat64List CSpectrumFluxCorrectionMeiksin::applyAdaptativeKernel(const TFloat64List& arr, 
                                                                 const Float64 z_center, 
                                                                 const std::shared_ptr<const CLSF>& lsf,
                                                                 const TFloat64List& lambdas)
{
    if(!arr.size()){
        throw GlobalException(INTERNAL_ERROR,"Cannot convolve: either kernel or array is empty. ");
    }

    Int32 n = arr.size(), Nhalf=-1;
    TFloat64List convolvedArr(arr);

    //determine the restframe convolution range, i.e., convolRange/(1+z_center)
    TFloat64Range convRange_rest(m_convolRange.GetBegin()/(1+z_center), m_convolRange.GetEnd()/(1+z_center));//conv range in restframe
    Int32 i_min = -1, i_max = -1; 
    bool ret = convRange_rest.getClosedIntervalIndices(lambdas, i_min, i_max, false); 
    if(!ret)
    {
        return convolvedArr;
    }
    Float64 tmp;
    for(Int32 i = i_min; i <= i_max; i++){
        Float64 lambda0 = lambdas[i];//lambda restframe
        //compute the adpative kernel at lambda0
        getLSFProfileVector(lambda0, z_center, lsf);//resulting kernel saved in m_kernel
        Nhalf = int(m_kernel.size()/2);
        if(!m_kernel.size()){
            throw GlobalException(INTERNAL_ERROR,"Cannot convolve: either kernel or array is empty. ");
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
/**
 * convolve only on m_convolRange/(1+zbin_meiksin), while keeping vector size
*/
void CSpectrumFluxCorrectionMeiksin::convolveAll(const std::shared_ptr<const CLSF>& lsf)
{
    //to be decided if we create two objects or we keep one m_correction
    m_corrections.resize(m_rawCorrections.size());

    TFloat64List meiksin_Bins = getSegmentsStartRedshiftList();
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
            TFloat64List a = applyAdaptativeKernel(m_rawCorrections[i].fluxcorr[j], z_center, lsf, m_corrections[i].lbda);
            m_corrections[i].fluxcorr.push_back(std::move(a));
        }
    }
    return;
}
