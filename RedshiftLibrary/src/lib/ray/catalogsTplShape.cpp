#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/linemodel/calibrationconfig.h>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/ray/catalogsTplShape.h>
#include <RedshiftLibrary/ray/linetags.h>

#include <algorithm>    // std::sort
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <string>
#include <fstream>
#include <iostream>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;
using namespace boost;


Bool CRayCatalogsTplShape::Init( std::string calibrationPath, std::string opt_tplratioCatRelPath, Int32 enableISMCalzetti)
{
    if(opt_tplratioCatRelPath.size()<1)
    {
      throw runtime_error("Unable to init the tpl-ratio catalog. Found empty relative path.");
    }

    bfs::path calibrationFolder( calibrationPath.c_str() );

    tplshapedcatalog_relpath = opt_tplratioCatRelPath;
    Log.LogInfo("    CatalogsTplShape - Loading tplshape catalog : %s", tplshapedcatalog_relpath.c_str());

    std::string dirPath = (calibrationFolder/tplshapedcatalog_relpath.c_str()).string();

    m_opt_dust_calzetti = enableISMCalzetti;
    //hardcoded fitting values for EBV, should be in the json
    Float64 ebmv_start=0.0;
    Float64 ebmv_step=0.1;
    Float64 ebmv_n=10;
    //Float64 ebmv_start=0.0;
    //Float64 ebmv_step=0.9;
    //Float64 ebmv_n=2;
    m_ismCorrectionCalzetti.Init(calibrationPath, ebmv_start, ebmv_step, ebmv_n);

    bool ret = Load(dirPath.c_str());
    if(!ret)
    {
        Log.LogError("    CatalogsTplShape: Unable to load the tpl-shape catalogs. aborting...");
        return false;
    }
    return true;
}


Bool CRayCatalogsTplShape::Load( const char* dirPath )
{
    // Clear current catalog list
    m_RayCatalogList.clear();
    m_RayCatalogNames.clear();
    m_ELvelocities.clear();
    m_ABSvelocities.clear();
    m_Priors.clear();
    m_IsmIndexes.clear();

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
    std::sort(tplshapeCatalogList.begin(),tplshapeCatalogList.end());
    Log.LogDebug("    CatalogsTplShape - Found %d tplshaped catalogs", tplshapeCatalogList.size());

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
        Log.LogError("    CatalogsTplShape - ERROR - unable to find velocities directory");
    }
    Log.LogInfo("    CatalogsTplShape - Found %d tplshaped velocities files", tplshapeVelocitiesList.size());

    //load the priors list for all the catalogs
    fs::path tplshapePriorsDir = tplshapeCatalogDir/"priors/";
    std::vector<std::string> tplshapePriorsList;
    if ( fs::exists(tplshapePriorsDir) && fs::is_directory(tplshapePriorsDir))
    {
      for( fs::directory_iterator dir_iter(tplshapePriorsDir) ; dir_iter != end_iter ; ++dir_iter)
      {
        if (fs::is_regular_file(dir_iter->status()) )
        {
          tplshapePriorsList.push_back(dir_iter->path().c_str());
        }
      }
    }else if ( fs::exists(tplshapePriorsDir) ) {
        Log.LogWarning("    CatalogsTplShape - ERROR - unable to find priors files");
    }
    Log.LogInfo("    CatalogsTplShape - Found %d tplshaped priors files", tplshapePriorsList.size());

    //load the priorsPz list for all the catalogs
    fs::path tplshapePriorsPzDir = tplshapeCatalogDir/"priors_pZ_Tplr/";
    std::vector<std::string> tplshapePriorsPzList;
    if ( fs::exists(tplshapePriorsPzDir) && fs::is_directory(tplshapePriorsPzDir))
    {
      for( fs::directory_iterator dir_iter(tplshapePriorsPzDir) ; dir_iter != end_iter ; ++dir_iter)
      {
        if (fs::is_regular_file(dir_iter->status()) )
        {
          tplshapePriorsPzList.push_back(dir_iter->path().c_str());
        }
      }
    }else if ( fs::exists(tplshapePriorsPzDir) ) {
        Log.LogWarning("    CatalogsTplShape - ERROR - unable to find priorsPz files");
    }
    Log.LogInfo("    CatalogsTplShape - Found %d tplshaped priorsPz files", tplshapePriorsPzList.size());


    //Load the linecatalog-tplshaped in the list
    Int32 successLoadPriors = true;
    Int32 nDustCoeffs=1;
    if(!m_opt_dust_calzetti)
    {
        nDustCoeffs = 1;
    }else{
        nDustCoeffs = m_ismCorrectionCalzetti.GetNPrecomputedDustCoeffs();
    }
    for(Int32 ktpl=0; ktpl<tplshapeCatalogList.size(); ktpl++)
    {
        //Loop on the EBMV dust coeff
        for(Int32 kDust=0; kDust<nDustCoeffs; kDust++)
        {
            CRayCatalog lineCatalog;
            try
            {
                lineCatalog.Load( tplshapeCatalogList[ktpl].c_str() );
            }
            catch (std::string& e)
            {
                Log.LogError("    CatalogsTplShape - Failed to load tplshape catalog: %s", e.c_str());
                throw runtime_error("    CatalogsTplShape - Failed to load tplshape catalog");
            }

            m_RayCatalogList.push_back(lineCatalog);
            m_IsmIndexes.push_back(kDust);

            fs::path name(tplshapeCatalogList[ktpl].c_str());
            m_RayCatalogNames.push_back(name.filename().c_str());

            //find the velocities-tplshaped corresponding to the tpl-shaped catalog: this is mandatory (continue/skip loading if velocities aren't found)
            Int32 kvel = -1;
            std::string tplname = name.filename().c_str();
            boost::replace_all( tplname, "_catalog.txt", "_velocities.txt");
            for(Int32 kv=0; kv<tplshapeVelocitiesList.size(); kv++)
            {
                std::string velname = tplshapeVelocitiesList[kv];
                std::size_t foundstra = velname.find(tplname.c_str());
                if (foundstra==std::string::npos){
                    continue;
                }
                kvel = kv;
            }

            if(kvel<0)
            {
                Log.LogError("    CatalogsTplShape - Failed to match tplshape-catalog with tplshape-velocities files: %s", tplname.c_str());
                return false;
            }
            m_ELvelocities.push_back(200.0);
            m_ABSvelocities.push_back(200.0);
            bool ret = LoadVelocities(tplshapeVelocitiesList[kvel].c_str(), m_ELvelocities.size()-1);
            if( !ret )
            {
                Log.LogError("    CatalogsTplShape - Failed to load tplshape velocities: %s", tplshapeVelocitiesList[kvel].c_str());
                return false;
            }

            //find the prior-tplshaped corresponding to the tpl-shaped catalog : this is optional (ignore if priors aren't found)
            m_Priors.push_back(1.0);
            Int32 kprior = -1;
            tplname = name.filename().c_str();
            boost::replace_all( tplname, "_catalog.txt", "_prior.txt");
            for(Int32 kp=0; kp<tplshapePriorsList.size(); kp++)
            {
                std::string priorname = tplshapePriorsList[kp];
                std::size_t foundstra = priorname.find(tplname.c_str());
                if (foundstra==std::string::npos){
                    continue;
                }
                kprior = kp;
            }

            if(kprior<0 || !successLoadPriors)
            {
                Log.LogDetail("    CatalogsTplShape - Failed to match tplshape-catalog with tplshape-prior files: %s", tplname.c_str());
                successLoadPriors=false;
            }
            if(successLoadPriors)
            {
                bool retPrior = LoadPrior(tplshapePriorsList[kprior].c_str(), m_Priors.size()-1);
                if( !retPrior )
                {
                    Log.LogError("    CatalogsTplShape - Failed to load tplshape prior: %s", tplshapePriorsList[kprior].c_str());
                    successLoadPriors=false;
                }else{
                    //make this prior tpl-ratio equiprobable with other EBMV values
                    if(nDustCoeffs>1)
                    {
                        m_Priors[m_Priors.size()-1] /= Float64(nDustCoeffs);
                    }
                }
            }
        }
    }
    if(!successLoadPriors) //if not all priors were successfully loaded, replace by cst prior
    {
        Float64 priorCST = 1./(Float64)(m_Priors.size());
        Log.LogDetail("    CatalogsTplShape - Failed to load tplshape prior, USING constant priors instead ! (p=%f)", priorCST);
        for(Int32 k=0; k<m_Priors.size(); k++)
        {
            for(Int32 kDust=0; kDust<nDustCoeffs; kDust++)
            {
                m_Priors[k] = priorCST;
            }
        }
    }
    Log.LogInfo("    CatalogsTplShape - Loaded %d tplshaped catalogs", m_RayCatalogList.size());

    return true;
}

bool CRayCatalogsTplShape::LoadVelocities( const char* filePath, Int32 k )
{
    Float64 elv=100.0;
    Float64 alv=300.0;

    std::ifstream file;
    file.open( filePath, std::ifstream::in );
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


    Log.LogDebug("    CatalogsTplShape - k=%d - Set elv=%.1f, alv=%.1f", k, elv, alv);
    m_ELvelocities[k] = elv;
    m_ABSvelocities[k] = alv;

    return true;
}

bool CRayCatalogsTplShape::LoadPrior( const char* filePath, Int32 k )
{
    Float64 prior=1.0;

    std::ifstream file;
    file.open( filePath, std::ifstream::in );
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
            prior = std::stod(line);
        }
        readNums++;
    }
    file.close();
    if(readNums!=1)
    {
        return false;
    }

    Log.LogDebug("    CatalogsTplShape - k=%d - Set prior=%f", k, prior);
    m_Priors[k] = prior;

    return true;
}

/**
 * @brief CRayCatalogsTplShape::GetRestLinesList
 * @param index
 * WARNING: ismCoeff not applied on the restlines provided by that function.
 */
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

const std::vector<Float64> & CRayCatalogsTplShape::getCatalogsPriors()
{
    return m_Priors;
}

std::string CRayCatalogsTplShape::GetCatalogName(Int32 idx)
{
    return m_RayCatalogNames[idx];
}

Float64 CRayCatalogsTplShape::GetIsmCoeff(Int32 idx)
{
    return m_ismCorrectionCalzetti.GetEbmvValue(m_IsmIndexes[idx]);
}

Int32 CRayCatalogsTplShape::GetIsmIndex(Int32 idx)
{
    return m_IsmIndexes[idx];
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
                Float64 restLambda = currentCatalogLineList[kL].GetPosition();
                Float64 dustCoeff = m_ismCorrectionCalzetti.getDustCoeff( m_IsmIndexes[iCatalog], restLambda);
                nominalAmp*=dustCoeff;
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

    //Now log the linesCorrespondingNominalAmp
    for( UInt32 iElts=0; iElts<LineModelElementList.m_Elements.size(); iElts++ )
    {
        for(Int32 k=0; k<GetCatalogsCount(); k++)
        {
            Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();
            for(UInt32 j=0; j<nRays; j++){
                Float64 ebv = m_ismCorrectionCalzetti.GetEbmvValue(m_IsmIndexes[k]);
                Float64 nomAmp = m_RayCatalogLinesCorrespondingNominalAmp[iElts][k][j];
                std::string lineName = LineModelElementList.m_RestRayList[LineModelElementList.m_Elements[iElts]->m_LineCatalogIndexes[j]].GetName();
                Log.LogDebug("    CatalogsTplShape - linesCorrespondingNominalAmp iElt=%d, iCatalog=%d, iLine=%d with name=%s, ebv=%f: NominalAmpFound = %e", iElts, k, j, lineName.c_str(), ebv, nomAmp);
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
        Float64 restLambda = currentCatalogLineList[kL].GetPosition();
        Float64 dustCoeff = m_ismCorrectionCalzetti.getDustCoeff( m_IsmIndexes[iCatalog], restLambda);
        nominalAmp*=dustCoeff;
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

Bool CRayCatalogsTplShape::SetLyaProfile(CLineModelElementList &LineModelElementList, Int32 iCatalog, bool forceLyaFitting)
{
    if(iCatalog<0){
        return false;
    }
    linetags ltags;
    std::string lyaTag = ltags.lya_em;
    //loop the amplitudes in the iLine_st catalog in order to find Lya
    CRayCatalog::TRayVector currentCatalogLineList = m_RayCatalogList[iCatalog].GetList();
    Int32 nLines = currentCatalogLineList.size();
    for(Int32 kL=0; kL<nLines; kL++)
    {
        if(! (currentCatalogLineList[kL].GetName()==lyaTag.c_str()))
        {
            continue;
        }
	CRay::TProfile targetProfile = currentCatalogLineList[kL].GetProfile();
    if(forceLyaFitting)
    {
        targetProfile = CRay::ASYMFIT;
    }

        //find line Lya in the elementList
        for( UInt32 iElts=0; iElts<LineModelElementList.m_Elements.size(); iElts++ )
        {
            //get the max nominal amplitude
            Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();
            for(UInt32 j=0; j<nRays; j++){

                if(LineModelElementList.m_RestRayList[LineModelElementList.m_Elements[iElts]->m_LineCatalogIndexes[j]].GetName() == lyaTag.c_str())
                {
                    TAsymParams asymParams = currentCatalogLineList[kL].GetAsymParams();
                    LineModelElementList.m_RestRayList[LineModelElementList.m_Elements[iElts]->m_LineCatalogIndexes[j]].SetProfile(targetProfile);
                    LineModelElementList.m_RestRayList[LineModelElementList.m_Elements[iElts]->m_LineCatalogIndexes[j]].SetAsymParams(asymParams);
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


