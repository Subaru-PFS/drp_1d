#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h>

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



namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;
using namespace boost;


CSpectrumFluxCorrectionCalzetti::CSpectrumFluxCorrectionCalzetti()
{
    m_LambdaMin = 100.0;
    m_LambdaMax = 99999.0;
}

CSpectrumFluxCorrectionCalzetti::~CSpectrumFluxCorrectionCalzetti()
{
    if(!calzettiInitFailed)
    {
        delete[] m_dataCalzetti;
        delete[] m_dataDustCoeff;
    }
}

Bool CSpectrumFluxCorrectionCalzetti::Init( std::string calibrationPath, Float64 ebmv_start, Float64 ebmv_step, Float64 ebmv_n)
{
    //load calzetti data
    bfs::path calibrationFolder( calibrationPath.c_str() );
    std::string filePath = (calibrationFolder/"ism"/"SB_calzetti.dl1.txt").string();
    std::ifstream file;
    file.open( filePath, std::ifstream::in );
    bool fileOpenFailed = file.rdstate() & std::ios_base::failbit;
    if(fileOpenFailed)
    {
        Log.LogError("ChisquareLog, unable to load the calzetti calib. file: %s... aborting!", filePath.c_str());
        calzettiInitFailed = true;
    }else
    {
        m_NdataCalzetti = 1e5;
        m_dataCalzetti = new Float64 [(int)m_NdataCalzetti]();

        std::string line;
        // Read file line by line
        Int32 kLine = 0;
        while( getline( file, line ) )
        {
            if( !boost::starts_with( line, "#" ) )
            {
                std::istringstream iss( line );
                Float64 x, y;
                iss >> x >> y;
                m_dataCalzetti[kLine] = y;
                kLine++;
                if(kLine>=m_NdataCalzetti)
                {
                    break;
                }
            }
        }
        file.close();

        //precomte the dust-coeff table
        m_nDustCoeff = ebmv_n;
        m_dustCoeffStep = ebmv_step;
        m_dustCoeffStart = ebmv_start;
        m_dataDustCoeff = new Float64[(int)(m_nDustCoeff*m_NdataCalzetti)]();

        for(Int32 kDust=0; kDust<m_nDustCoeff; kDust++)
        {

            Float64 coeffEBMV = GetEbmvValue(kDust);
            for(Int32 kCalzetti=0; kCalzetti<m_NdataCalzetti; kCalzetti++)
            {
                m_dataDustCoeff[Int32(kDust*m_NdataCalzetti+kCalzetti)] = pow(10.0, -0.4*m_dataCalzetti[kCalzetti]*coeffEBMV);
            }

        }

        calzettiInitFailed = false;
    }
    return true;
}

Float64 CSpectrumFluxCorrectionCalzetti::GetEbmvValue(Int32 k)
{
    Float64 coeffEBMV = m_dustCoeffStart + m_dustCoeffStep*(Float64)k;
    return coeffEBMV;
}

Float64 CSpectrumFluxCorrectionCalzetti::getDustCoeff( Int32 kDust, Float64 restLambda )
{
    Float64 coeffDust = 1.0;
    if(restLambda >= m_LambdaMin && restLambda < m_LambdaMax)
    {
        Int32 kCalzetti = Int32(restLambda-100.0);
        coeffDust = m_dataDustCoeff[Int32(kDust*m_NdataCalzetti+kCalzetti)];
    }
    return coeffDust;
}

/* @brief CSpectrumFluxCorrectionCalzetti::getDustCoeff: get the dust coeff at a fixed resolution of 1A
* @param dustCoeff
* @param maxLambda
* @return
*/
const Float64*  CSpectrumFluxCorrectionCalzetti::getDustCoeff(Float64 dustCoeff, Float64 maxLambda)
{
    //find kDust
    Int32 idxDust = -1;
    for(Int32 kDust=0; kDust<m_nDustCoeff; kDust++)
    {
        Float64 coeffEBMV =  GetEbmvValue(kDust);
        if(dustCoeff==coeffEBMV)
        {
            idxDust = kDust;
            break;
        }
    }
    if(idxDust<0)
    {
        return 0;
    }

    Int32 nSamples = maxLambda+1; //+1 for security
    Float64* dustCoeffs = new Float64 [(int)nSamples]();


    for(Int32 kl=0; kl<nSamples; kl++)
    {
        Float64 restLambda = kl;
        Float64 coeffDust = 1.0;
        if(restLambda >= 100.0)
        {
            Int32 kCalzetti = Int32(restLambda-100.0);
            coeffDust = m_dataDustCoeff[Int32(idxDust*m_NdataCalzetti+kCalzetti)];
        }
        dustCoeffs[kl] = coeffDust;
    }
    return dustCoeffs;
}

Int32 CSpectrumFluxCorrectionCalzetti::GetNPrecomputedDustCoeffs()
{
    return m_nDustCoeff;
}

Float64 CSpectrumFluxCorrectionCalzetti::GetLambdaMin()
{
    return m_LambdaMin;
}

Float64 CSpectrumFluxCorrectionCalzetti::GetLambdaMax()
{
    return m_LambdaMax;
}


