#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/catalogsTplShape.h>

#include <algorithm>    // std::sort
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <string>
#include <fstream>
#include <iostream>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;
using namespace boost;


CRayCatalogsTplShape::CRayCatalogsTplShape()
{
    // old relpath = "linecatalogs_tplshape_ExtendedTemplatesMarch2016_v2_20160916_B10I2_mod"
    //tplshapedcatalog_relpath = "linecatalogs_tplshape_ExtendedTemplatesMarch2016_B13B_mod20170110";
    //tplshapedcatalog_relpath = "linecatalogs_tplshape_ExtendedTemplatesMarch2016_B13D_mod2";

    bfs::path tplshapeRelPath( "linecatalogs_tplshapes" );
    tplshapedcatalog_relpath = (tplshapeRelPath/"linecatalogs_tplshape_ExtendedTemplatesJan2017v3_20170524_B13F_v1").string();
    Log.LogInfo( "CRayCatalogsTplShape - Loaded tplshape catalog : %s", tplshapedcatalog_relpath.c_str());
}

CRayCatalogsTplShape::~CRayCatalogsTplShape()
{

}

Bool CRayCatalogsTplShape::SetTplctlgRelPath( const char* relPath )
{
    tplshapedcatalog_relpath = relPath;
    return true;
}

Bool CRayCatalogsTplShape::Init( std::string calibrationPath)
{
    bfs::path calibrationFolder( calibrationPath.c_str() );
    //std::string dirPath = (calibrationFolder.append( tplshapedcatalog_relpath.c_str() )).string();
    std::string dirPath = (calibrationFolder/tplshapedcatalog_relpath.c_str()).string();

    bool ret = Load(dirPath.c_str());
    if(!ret)
    {
        Log.LogError("Unable to load the tpl-shape catalogs. aborting...");
        return false;
    }
    return true;
}


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
    if(tplshapeCatalogList.size()<1)
    {
        return false;
    }
    Log.LogDebug( "CRayCatalogsTplShape - Found %d tplshaped catalogs", tplshapeCatalogList.size());

    //load the velocities list for all the catalogs
    fs::path tplshapeVelocitiesDir = tplshapeCatalogDir/"velocities/";
    std::vector<std::string> tplshapeVelocitiesList;
    if ( fs::exists(tplshapeVelocitiesDir) && fs::is_directory(tplshapeVelocitiesDir))
    {
      for( fs::directory_iterator dir_iter(tplshapeVelocitiesDir) ; dir_iter != end_iter ; ++dir_iter)
      {
        if (fs::is_regular_file(dir_iter->status()) )
        {
          tplshapeVelocitiesList.push_back(dir_iter->path().c_str());
        }
      }
    }else{
        Log.LogError( "CRayCatalogsTplShape - ERROR - unable to find velocities directory");
    }
    Log.LogInfo( "CRayCatalogsTplShape - Found %d tplshaped velocities files", tplshapeVelocitiesList.size());



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

        fs::path name(tplshapeCatalogList[k].c_str());
        m_RayCatalogNames.push_back(name.filename().c_str());

        //find the velocities-tplshaped corresponding to the tpl-shaped catalog
        Int32 kvel = -1;
        std::string tplname = name.filename().c_str();
        boost::replace_all( tplname, "_catalog.txt", "_velocities.txt");
        for(Int32 k=0; k<tplshapeVelocitiesList.size(); k++)
        {
            std::string velname = tplshapeVelocitiesList[k];
            std::size_t foundstra = velname.find(tplname.c_str());
            if (foundstra==std::string::npos){
                continue;
            }
            kvel = k;
        }

        if(kvel<0)
        {
            Log.LogError( "Failed to match tplshape-catalog with tplshape-velocities files: %s", tplname.c_str());
            return false;
        }
        m_ELvelocities.push_back(200.0);
        m_ABSvelocities.push_back(200.0);
        bool ret = LoadVelocities(tplshapeVelocitiesList[kvel].c_str(), k);
        if( !ret )
        {
            Log.LogError( "Failed to load tplshape velocities: %s", tplshapeVelocitiesList[kvel].c_str());
            return false;
        }


    }
    Log.LogInfo( "CRayCatalogsTplShape - Loaded %d tplshaped catalogs", m_RayCatalogList.size());

    return true;
}

bool CRayCatalogsTplShape::LoadVelocities( const char* filePath, Int32 k )
{
    Float64 elv=100.0;
    Float64 alv=300.0;

    ifstream file;
    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit ){
        return false;
    }
    string line;

    // Read file line by line
    Int32 readNums = 0;
    while( getline( file, line ) )
    {
        if(readNums==0)
        {
            elv = std::stod(line);
        }else if(readNums==1)
        {
            alv = std::stod(line);
        }
        readNums++;
    }
    file.close();
    if(readNums!=2)
    {
        return false;
    }


    Log.LogDebug( "CRayCatalogsTplShape k=%d - Set elv=%.1f, alv=%.1f", k, elv, alv);
    m_ELvelocities[k] = elv;
    m_ABSvelocities[k] = alv;

    return true;
}

CRayCatalog::TRayVector CRayCatalogsTplShape::GetRestLinesList( const Int32 index )
{
    Int32 typeFilter=-1;
    Int32 forceFilter=-1;

    CRayCatalog::TRayVector restRayList = m_RayCatalogList[index].GetFilteredList( typeFilter, forceFilter);
    return restRayList;
}

Int32 CRayCatalogsTplShape::GetCatalogsCount()
{
    return m_RayCatalogList.size();
}

std::string CRayCatalogsTplShape::GetCatalogName(Int32 idx)
{
    return m_RayCatalogNames[idx];
}

Bool CRayCatalogsTplShape::GetCatalogVelocities(Int32 idx, Float64& elv, Float64& alv )
{
    elv = m_ELvelocities[idx];
    alv = m_ABSvelocities[idx];
    return true;
}


Bool CRayCatalogsTplShape::InitLineCorrespondingAmplitudes(CLineModelElementList &LineModelElementList)
{
    //first set all corresponding amplitudes to 0.0;
    for( UInt32 iElts=0; iElts<LineModelElementList.m_Elements.size(); iElts++ )
    {
        std::vector<std::vector<Float64>> thisElementLinesCorresp;

        Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();

        for(Int32 k=0; k<GetCatalogsCount(); k++)
        {
            std::vector<Float64> thisCatLinesCorresp;
            for(UInt32 j=0; j<nRays; j++){
                thisCatLinesCorresp.push_back(0.0);
            }
            thisElementLinesCorresp.push_back(thisCatLinesCorresp);
        }
        m_RayCatalogLinesCorrespondingNominalAmp.push_back(thisElementLinesCorresp);

        //now set the non-zero amp correspondences
        for(Int32 iCatalog=0; iCatalog<GetCatalogsCount(); iCatalog++)
        {
            CRayCatalog::TRayVector currentCatalogLineList = m_RayCatalogList[iCatalog].GetList();
            for(Int32 kL=0; kL<currentCatalogLineList.size(); kL++)
            {
                Float64 nominalAmp = currentCatalogLineList[kL].GetNominalAmplitude();
                //find line in the elementList
                Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();
                for(UInt32 j=0; j<nRays; j++){

                    if(LineModelElementList.m_RestRayList[LineModelElementList.m_Elements[iElts]->m_LineCatalogIndexes[j]].GetName() == currentCatalogLineList[kL].GetName())
                    {
                        m_RayCatalogLinesCorrespondingNominalAmp[iElts][iCatalog][j]=nominalAmp;
                    }
                }
            }
        }
    }

    return 0;
}

/**
 * @brief CRayCatalogsTplShape::SetMultilineNominalAmplitudesFast
 * This method sets the linemodel unique elt nominal amplitudes to the corresponding value of the iCatalog st catalog.
 * INFO: fast method, InitLineCorrespondence() should have been called previously with the same LineModelElementList arg.
 * @param LineModelElementList
 * @param iCatalog
 * @return
 */
Bool CRayCatalogsTplShape::SetMultilineNominalAmplitudesFast(CLineModelElementList &LineModelElementList, Int32 iCatalog)
{
    if(iCatalog<0){
        return false;
    }
    Float64 nominalAmp = 0.0;
    for( UInt32 iElts=0; iElts<LineModelElementList.m_Elements.size(); iElts++ )
    {
        Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();
        for(UInt32 j=0; j<nRays; j++){
            nominalAmp = m_RayCatalogLinesCorrespondingNominalAmp[iElts][iCatalog][j];
            LineModelElementList.m_Elements[iElts]->SetNominalAmplitude(j, nominalAmp);
        }
    }
    return true;
}

/**
 * @brief CRayCatalogsTplShape::SetMultilineNominalAmplitudes
 * This method sets the linemodel unique elt nominal amplitudes to the corresponding value of the iCatalog st catalog.
 * INFO: slow method
 * @param LineModelElementList
 * @param iCatalog
 * @return
 */
Bool CRayCatalogsTplShape::SetMultilineNominalAmplitudes(CLineModelElementList &LineModelElementList, Int32 iCatalog)
{
    //first set all amplitudes to 0.0
    for( UInt32 iElts=0; iElts<LineModelElementList.m_Elements.size(); iElts++ )
    {
        //get the max nominal amplitude
        Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();
        for(UInt32 j=0; j<nRays; j++){
            LineModelElementList.m_Elements[iElts]->SetNominalAmplitude(j, 0.0);
        }
    }

    //loop the amplitudes in the iLine_st catalog
    CRayCatalog::TRayVector currentCatalogLineList = m_RayCatalogList[iCatalog].GetList();
    Int32 nLines = currentCatalogLineList.size();
    for(Int32 kL=0; kL<nLines; kL++)
    {
        Float64 nominalAmp = currentCatalogLineList[kL].GetNominalAmplitude();
        //find line in the elementList
        for( UInt32 iElts=0; iElts<LineModelElementList.m_Elements.size(); iElts++ )
        {
            //get the max nominal amplitude
            Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();
            for(UInt32 j=0; j<nRays; j++){

                if(LineModelElementList.m_RestRayList[LineModelElementList.m_Elements[iElts]->m_LineCatalogIndexes[j]].GetName() == currentCatalogLineList[kL].GetName())
                {
                    LineModelElementList.m_Elements[iElts]->SetNominalAmplitude(j, nominalAmp);
                }

            }


        }

    }
    return true;
}

Bool CRayCatalogsTplShape::SetLyaProfile(CLineModelElementList &LineModelElementList, Int32 iCatalog)
{
    if(iCatalog<0){
        return false;
    }
    std::string lyaTag = "LyAE";
    //loop the amplitudes in the iLine_st catalog in order to find Lya
    CRayCatalog::TRayVector currentCatalogLineList = m_RayCatalogList[iCatalog].GetList();
    Int32 nLines = currentCatalogLineList.size();
    for(Int32 kL=0; kL<nLines; kL++)
    {
        if(! (currentCatalogLineList[kL].GetName()==lyaTag.c_str()))
        {
            continue;
        }
        std::string targetProfile = currentCatalogLineList[kL].GetProfile();

        //find line Lya in the elementList
        for( UInt32 iElts=0; iElts<LineModelElementList.m_Elements.size(); iElts++ )
        {
            //get the max nominal amplitude
            Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();
            for(UInt32 j=0; j<nRays; j++){

                if(LineModelElementList.m_RestRayList[LineModelElementList.m_Elements[iElts]->m_LineCatalogIndexes[j]].GetName() == lyaTag.c_str())
                {

                    LineModelElementList.m_RestRayList[LineModelElementList.m_Elements[iElts]->m_LineCatalogIndexes[j]].SetProfile(targetProfile);
                    break;
                }

            }


        }

    }
    return true;
}


/**
 * \brief Calculates the best fit between the linemodel fitted amplitudes and the tplShaped catalogs: (for lm-rigidity=tplcorr)
 *
 **/
Float64 CRayCatalogsTplShape::GetBestFit( const CRayCatalog::TRayVector& restRayList, std::vector<Float64> fittedAmplitudes, std::vector<Float64> fittedErrors, std::vector<Float64>& amplitudesCorrected, std::string& bestTplName  )
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
            if(fittedAmplitudes[iRestRay]>=0 && fittedErrors[iRestRay]>0)
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
                bestTplName = m_RayCatalogNames[iCatalogs];
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
            iTplAmps++;
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
    Float64 sumCross= 0.0;
    Float64 sumTPL2 = 0.0;

    for(Int32 k=0; k<N; k++)
    {
        Float64 err2 = 1.0 / (errLM[k] * errLM[k]);

        sumCross+=ampsLM[k]*ampsTPL[k]*err2;
        sumTPL2+=ampsTPL[k]*ampsTPL[k]*err2;

        sumLM += ampsLM[k]*err2;
        sumTPL += ampsTPL[k]*err2;
    }
//    if ( sumLM==0 || sumTPL==0 )
//    {
//        return -1.0;
//    }
//    Float64 ampl = sumLM / sumTPL;
    if ( sumCross==0 || sumTPL2==0 )
    {
        return -1.0;
    }
    Float64 ampl = sumCross / sumTPL2;

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


