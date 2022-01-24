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
#include "RedshiftLibrary/linemodel/calibrationconfig.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/ray/catalogsTplShape.h"
#include "RedshiftLibrary/ray/linetags.h"

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


bool CRayCatalogsTplShape::Init( std::string calibrationPath, 
                                std::string opt_tplratioCatRelPath, 
                                Int32 enableISMCalzetti, 
                                std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti,
                                Float64 nsigmasupport)
{
    m_nsigmasupport = nsigmasupport;
    m_opt_dust_calzetti = enableISMCalzetti;
    m_ismCorrectionCalzetti = ismCorrectionCalzetti;

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

    CRayCatalog::TRayVector restRayList = m_lineRatioCatalogs[index].GetFilteredList( typeFilter, forceFilter);
    return restRayList;
}

Int32 CRayCatalogsTplShape::GetCatalogsCount()
{
    return m_lineRatioCatalogs.size();
}

std::vector<Float64>  CRayCatalogsTplShape::getCatalogsPriors()
{
  std::vector<Float64> ret;
  for(CLineRatioCatalog cat : m_lineRatioCatalogs) ret.push_back(cat.getPrior());
  return ret;
}

std::string CRayCatalogsTplShape::GetCatalogName(Int32 idx)
{
  return m_lineRatioCatalogs[idx].getName();
}

Float64 CRayCatalogsTplShape::GetIsmCoeff(Int32 idx)
{
  return m_ismCorrectionCalzetti->GetEbmvValue(GetIsmIndex(idx));
}

Int32 CRayCatalogsTplShape::GetIsmIndex(Int32 idx)
{
  return m_lineRatioCatalogs[idx].getIsmIndex();
}

bool CRayCatalogsTplShape::GetCatalogVelocities(Int32 idx, Float64& elv, Float64& alv )
{
  //TODO generic velocity groups : there should not be hardcoded values, this should return a map
  elv = m_lineRatioCatalogs[idx].getVelocity("em_vel");
    alv = m_lineRatioCatalogs[idx].getVelocity("abs_vel");
    return true;
}


bool CRayCatalogsTplShape::InitLineCorrespondingAmplitudes(const CLineModelElementList &LineModelElementList)
{
    //first set all corresponding amplitudes to 0.0;
    for( UInt32 iElts=0; iElts<LineModelElementList.size(); iElts++ )
    {
        Int32 nRays = LineModelElementList[iElts]->GetSize();
        TFloat64List thisCatLinesCorresp(nRays, 0.0); //is nRays cte among Elts?
        std::vector<TFloat64List> thisElementLinesCorresp(GetCatalogsCount(), thisCatLinesCorresp);

        m_RayCatalogLinesCorrespondingNominalAmp.push_back(thisElementLinesCorresp);

        //now set the non-zero amp correspondences
        for(Int32 iCatalog=0; iCatalog<GetCatalogsCount(); iCatalog++)
        {
            CRayCatalog::TRayVector currentCatalogLineList = m_lineRatioCatalogs[iCatalog].GetList();
            for(Int32 kL=0; kL<currentCatalogLineList.size(); kL++)
            {
                Float64 nominalAmp = currentCatalogLineList[kL].GetNominalAmplitude();
                Float64 restLambda = currentCatalogLineList[kL].GetPosition();
                Float64 dustCoeff = m_ismCorrectionCalzetti->GetDustCoeff( GetIsmIndex(iCatalog), restLambda);
                nominalAmp*=dustCoeff;
                //find line in the elementList
                Int32 nRays = LineModelElementList[iElts]->GetSize();
                for(UInt32 j=0; j<nRays; j++){

                    if(LineModelElementList[iElts]->m_Rays[j].GetName() == currentCatalogLineList[kL].GetName())
                    {
                        m_RayCatalogLinesCorrespondingNominalAmp[iElts][iCatalog][j]=nominalAmp;
                    }
                }
            }
        }
    }

    //Now log the linesCorrespondingNominalAmp
    for( UInt32 iElts=0; iElts<LineModelElementList.size(); iElts++ )
    {
        for(Int32 k=0; k<GetCatalogsCount(); k++)
        {
            Int32 nRays = LineModelElementList[iElts]->GetSize();
            for(UInt32 j=0; j<nRays; j++){
	      Float64 ebv = m_ismCorrectionCalzetti->GetEbmvValue(GetIsmIndex(k));
                Float64 nomAmp = m_RayCatalogLinesCorrespondingNominalAmp[iElts][k][j];
                std::string lineName = LineModelElementList[iElts]->m_Rays[j].GetName();
                Log.LogDebug("    CatalogsTplShape - linesCorrespondingNominalAmp iElt=%d, iCatalog=%d, iLine=%d with name=%s, ebv=%f: NominalAmpFound = %e", iElts, k, j, lineName.c_str(), ebv, nomAmp);
            }
        }
    }

    return 0;
}

const CRayCatalog& CRayCatalogsTplShape::GetCatalog(Int32 iCatalog)
{
    return m_lineRatioCatalogs[iCatalog];
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
    for(UInt32 iCatalogs=0; iCatalogs<m_lineRatioCatalogs.size(); iCatalogs++)
    {
        CRayCatalog::TRayVector currentCatalogLineList = m_lineRatioCatalogs[iCatalogs].GetList();

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
            if(fit>0.0 && !std::isnan(fit) && (fit<coeffMin || coeffMin==-1))
            {
                coeffMin = fit;
                bestFitAmplitudes = ampsCorrected;
                bestTplName = GetCatalogName(iCatalogs);
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


