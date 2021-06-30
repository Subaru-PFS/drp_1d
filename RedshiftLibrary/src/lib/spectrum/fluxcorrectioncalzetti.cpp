#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

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
        std::string line;
        // Read file line by line
        while( getline( file, line ) )
        {
            if( !boost::starts_with( line, "#" ) )
            {
                std::istringstream iss( line );
                Float64 x, y;
                iss >> x >> y;
                m_dataCalzetti.push_back(y);
            }
        }
        file.close();
        
        //precomte the dust-coeff table
        m_nEbmvCoeff = ebmv_n;
        m_EbmvCoeffStep = ebmv_step;
        m_EbmvCoeffStart = ebmv_start;
        m_dataDustCoeff.resize(m_nEbmvCoeff*m_dataCalzetti.size());

        for(Int32 kDust=0; kDust<m_nEbmvCoeff; kDust++)
        {

            Float64 coeffEBMV = GetEbmvValue(kDust);
            for(Int32 kCalzetti=0; kCalzetti<m_dataCalzetti.size(); kCalzetti++)
            {
                m_dataDustCoeff[kDust*m_dataCalzetti.size()+kCalzetti] = pow(10.0, -0.4*m_dataCalzetti[kCalzetti]*coeffEBMV);
            }

        }

        calzettiInitFailed = false;
    }
    return true;
}

Float64 CSpectrumFluxCorrectionCalzetti::GetEbmvValue(Int32 k) const
{
    Float64 coeffEBMV = m_EbmvCoeffStart + m_EbmvCoeffStep*(Float64)k;
    return coeffEBMV;
}

Int32 CSpectrumFluxCorrectionCalzetti::GetEbmvIndex(Float64 value) const
{
    Int32 kEbmv = round((value - m_EbmvCoeffStart)/m_EbmvCoeffStep);
    return kEbmv;
}

Float64 CSpectrumFluxCorrectionCalzetti::GetDustCoeff( Int32 kDust, Float64 restLambda ) const 
{
    Float64 coeffDust = 1.0;
    if(restLambda >= m_LambdaMin && restLambda < m_LambdaMax)
    {
        Int32 kCalzetti = Int32(restLambda-100.0);
        coeffDust = m_dataDustCoeff[kDust*m_dataCalzetti.size()+kCalzetti];
    }
    return coeffDust;
}

Int32 CSpectrumFluxCorrectionCalzetti::GetNPrecomputedEbmvCoeffs() const
{
    return m_nEbmvCoeff;
}

Float64 CSpectrumFluxCorrectionCalzetti::GetLambdaMin() const
{
    return m_LambdaMin;
}

Float64 CSpectrumFluxCorrectionCalzetti::GetLambdaMax() const
{
    return m_LambdaMax;
}


