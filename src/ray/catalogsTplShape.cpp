#include <epic/core/log/log.h>
#include <epic/redshift/ray/catalogsTplShape.h>

#include <algorithm>    // std::sort
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <fstream>
#include <iostream>

using namespace NSEpic;
using namespace std;
using namespace boost;


CRayCatalogsTplShape::CRayCatalogsTplShape()
{

}

CRayCatalogsTplShape::~CRayCatalogsTplShape()
{

}

Bool CRayCatalogsTplShape::Init()
{
    std::string dirPath = "/home/aschmitt/data/vuds/VUDS_flag3_4/amazed/linecatalogs/linecatalogs_tplshape_ExtendedTemplatesMarch2016_v2_20160730_B10H";
    Load(dirPath.c_str());
    //Bool catalogsAreAligned = AreCatalogsAligned(restRayList, typeFilter, forceFilter);
    return true;
}

//**
// * \brief DEPRECATED
// * loops on the catalogstplshape and check if every one is aligned with the variable restRayList
// * NB: aligned means that the same line is accessed with a given list index
// * Non alignment can happen if the amazed catalog loaded (restRayList) is different than the tplShapedCatalogs.
// **/
//Bool CRayCatalogsTplShape::AreCatalogsAligned( const CRayCatalog::TRayVector& restRayList, Int32 typeFilter, Int32 forceFilter )
//{
//    Bool aligned=true;
//    for(UInt32 iCatalogs=0; iCatalogs<m_RayCatalogList.size(); iCatalogs++)
//    {
//        CRayCatalog::TRayVector currentCatalogLineList = m_RayCatalogList[iCatalogs].GetFilteredList( typeFilter, forceFilter);

//        for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
//        {
//            if(restRayList[iRestRay].GetType() != currentCatalogLineList[iRestRay].GetType()){
//                aligned=false;
//                break;
//            }
//            if(restRayList[iRestRay].GetName() != currentCatalogLineList[iRestRay].GetName()){
//                aligned=false;
//                break;
//            }
//            Float64 thres=1e-4;
//            Float64 positionDiff = abs(restRayList[iRestRay].GetPosition()-currentCatalogLineList[iRestRay].GetPosition());
//            if( positionDiff>thres ){
//                aligned=false;
//                break;
//            }
//        }
//        if(!aligned)
//        {
//            break;
//        }
//    }
//    return aligned;
//}


Bool CRayCatalogsTplShape::Load( const char* dirPath )
{

    // Clear current catalog list
    m_RayCatalogList.clear();

    //load the catalogs list from the files in the tplshape-catalogs folder : tplshapeCatalogDir
    namespace fs = boost::filesystem;
    fs::path tplshapeCatalogDir(dirPath);

    fs::directory_iterator end_iter;
    std::vector<std::string> tplshapeCatalogList;
    if ( fs::exists(tplshapeCatalogDir) && fs::is_directory(tplshapeCatalogDir))
    {
      for( fs::directory_iterator dir_iter(tplshapeCatalogDir) ; dir_iter != end_iter ; ++dir_iter)
      {
        if (fs::is_regular_file(dir_iter->status()) )
        {
          tplshapeCatalogList.push_back(dir_iter->path().c_str());
        }
      }
    }
    Log.LogDebug( "CRayCatalogsTplShape - Found %d tplshaped catalogs", tplshapeCatalogList.size());


    //Load the linecatalog-tplshaped in the list
    for(Int32 k=0; k<tplshapeCatalogList.size(); k++)
    {
        CRayCatalog lineCatalog;
        Bool rValue = lineCatalog.Load( tplshapeCatalogList[k].c_str() );
        if( !rValue )
        {
            Log.LogError( "Failed to load tplshape catalog: %s", tplshapeCatalogList[k].c_str());
            continue;
        }
        m_RayCatalogList.push_back(lineCatalog);
    }
    Log.LogDebug( "CRayCatalogsTplShape - Loaded %d tplshaped catalogs", m_RayCatalogList.size());

    return true;
}

/**
 * \brief Calculates the best fit between the linemodel fitted amplitudes and the tplShaped catalogs
 *
 **/
Float64 CRayCatalogsTplShape::GetBestFit( const CRayCatalog::TRayVector& restRayList, std::vector<Float64> fittedAmplitudes, std::vector<Float64> fittedErrors, std::vector<Float64>& amplitudesCorrected  )
{
    Float64 coeffMin = -1;
    std::vector<Int32> mask;
    std::vector<Float64> bestFitAmplitudes;
    //bestFitAmplitudes.resize(restRayList.size());
    for(UInt32 iCatalogs=0; iCatalogs<m_RayCatalogList.size(); iCatalogs++)
    {
        CRayCatalog::TRayVector currentCatalogLineList = m_RayCatalogList[iCatalogs].GetList();

        //create the amplitude float vectors
        std::vector<Float64> tplshapeAmplitudes;
        std::vector<Float64> linemodelAmplitudes;
        std::vector<Float64> linemodelErrors;

        mask.resize(fittedAmplitudes.size());

        for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
        {
            if(fittedAmplitudes[iRestRay]>=0)
            {
                mask[iRestRay]=1;
                Int32 iTplshapeRayFound = -1;
                for(UInt32 itplshapeRay=0; itplshapeRay<currentCatalogLineList.size(); itplshapeRay++)
                {
                    std::string tplshapeRayName =  currentCatalogLineList[itplshapeRay].GetName();
                    std::string restRayName = restRayList[iRestRay].GetName();
                    if ( restRayName==tplshapeRayName ){
                        iTplshapeRayFound = itplshapeRay;
                        break;
                    }
                }

                if ( iTplshapeRayFound < 0 )
                {
                    tplshapeAmplitudes.push_back(0.0);
                }else
                {
                    Float64 amp = currentCatalogLineList[iTplshapeRayFound].GetNominalAmplitude();
                    tplshapeAmplitudes.push_back(amp);
                }
                linemodelAmplitudes.push_back(fittedAmplitudes[iRestRay]);
                linemodelErrors.push_back(fittedErrors[iRestRay]);
            }else
            {
                mask[iRestRay]=0;
            }
        }

        if(linemodelAmplitudes.size()>1 && linemodelAmplitudes.size()==tplshapeAmplitudes.size())
        {
            std::vector<Float64> ampsCorrected;
            ampsCorrected.resize(linemodelAmplitudes.size());
            Float64 fit = GetFit(linemodelAmplitudes, linemodelErrors, tplshapeAmplitudes, ampsCorrected);
            if(fit>0.0 && !boost::math::isnan(fit) && (fit<coeffMin || coeffMin==-1))
            {
                coeffMin = fit;
                bestFitAmplitudes = ampsCorrected;
            }
        }

    }
    //coeff min normalization
    //coeffMin = sqrt(coeffMin);
    //coeffMin /= 1.0;

    //fill the corrected amplitudes vector
    Int32 iTplAmps = 0;
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        if(mask[iRestRay]>0 && coeffMin>=0)
        {
            amplitudesCorrected[iRestRay]=bestFitAmplitudes[iTplAmps];
        }else
        {
            amplitudesCorrected[iRestRay]=-1.0;
        }

    }

    return coeffMin;
}

Float64 CRayCatalogsTplShape::GetFit( std::vector<Float64> ampsLM, std::vector<Float64> errLM, std::vector<Float64> ampsTPL, std::vector<Float64>& ampsCorrected )
{

    Float64 N = ampsLM.size();
//    // Normalize AmpsLM first
//    Float64 normalizeCoeff = 0.0;
//    for(Int32 k=0; k<N; k++)
//    {
//        normalizeCoeff+= ampsLM[k];
//    }
//    for(Int32 k=0; k<N; k++)
//    {
//        ampsLM[k] = ampsLM[k]/normalizeCoeff;
//    }
    Float64 normalizeCoeff = 1.0;

    //estimate fitting amplitude
    Float64 sumLM = 0.0;
    Float64 sumTPL = 0.0;
    for(Int32 k=0; k<N; k++)
    {
        Float64 err2 = 1.0 / (errLM[k] * errLM[k]);
        sumLM += ampsLM[k]*err2;
        sumTPL += ampsTPL[k]*err2;
    }
    if ( sumLM==0 || sumTPL==0 )
    {
        return -1.0;
    }
    Float64 ampl = sumLM / sumTPL;

    Float64 fit=0.0;
    Float64 diff;
    for(Int32 k=0; k<N; k++)
    {
        Float64 err2 = 1.0 / (errLM[k] * errLM[k]);
        diff = ampsLM[k]-ampl*ampsTPL[k];
        fit += diff*diff*err2;
    }

    //fill the amps_corrected vector
    for(Int32 k=0; k<N; k++)
    {
        ampsCorrected[k] = ampsTPL[k]*ampl*normalizeCoeff;
    }

    return fit;
}


